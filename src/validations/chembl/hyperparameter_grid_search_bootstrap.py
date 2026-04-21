# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 10:40:00 2026

Bootstrap hyperparameter search for ChEMBL IC50 vs connectivity score validation.

This script is intentionally built on top of final_step_correlations_IC50_CS.py.
It reuses its data-loading, collapsing, classification, and plotting helpers, then:

1. evaluates MANY fixed hyperparameter combinations under stratified bootstrap resampling
2. stores bootstrap replicate metrics per configuration
3. summarizes precision / recall / F1 across configurations
4. plots a matrix of fig3-like scatter panels, one panel per configuration
5. reports also the metric on the full dataset for each configuration

Important scope note
--------------------
This script only tunes hyperparameters that can be evaluated from already-available
CS / IC50 tables. Hyperparameters that require recomputing the upstream CS itself
(e.g. LM vs all genes, DEG vs MITH, pathway-level vs gene-level signatures, different
CS runs) can still be included, but ONLY if the corresponding CS files / run IDs already
exist and are passed through the grid.

Recommended use
---------------
Use bootstrap here as a resampling-based estimate of metric variability across drugs.
This script does not train a model; it evaluates fixed configurations.

@Author: L-F-S
"""

from __future__ import annotations

import itertools
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2] if len(HERE.parents) >= 3 else HERE.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from conf import CS_DIR, CS_OUT, DATA_DIR, DISEASE,IC50_DRUG_COLLAPSE_METHOD,IC50_EFF_TH,\
    IC50_ONLY, IMG_DIR, LINCS_METADATA_PATH,LOGS_DIR, cell_line,cell_line_run_name,\
    cs_filename,cs_log_filename,cs_mith,cs_on_LM,disease_run_name,\
    ic50_file, selected_cs_run_id
from logger import append_run_metadata
from loader import load_IC50, load_drug_rankings
from validations.chembl.final_step_correlations_IC50_CS import add_pert_iname_to_cs, classify_ic50_vs_cs,\
    collapse_profiles_to_drug, plot_binchen_fig3_style, resolve_cs_run_id


# -----------------------------------------------------------------------------
# config object
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class EvalConfig:
    """
    Container for one hyperparameter configuration.
    """
    name: str
    disease: str
    cell_line: str
    cs_out_dir: Path
    cs_run_id: str
    cs_filename: str

    cell_line_IC50: str
    cs_drug_colname: str = "pert_iname"
    ic50_drug_colname: str = "pert_iname"
    cs_score_col: str = "connectivity_score"
    ic50_score_col: str = "standard_value"
    cs_drug_collapse_method: str = "best"
    ic50_drug_collapse_method: str = "median"
    ic50_only: bool = True

    cs_threshold: float = -1.5
    ic50_threshold: float = 10.0

    lincs_metadata_path: Path = LINCS_METADATA_PATH
    ic50_file: Path = ic50_file


# -----------------------------------------------------------------------------
# utilities
# -----------------------------------------------------------------------------

def expand_grid(param_grid):
    """
    Create cartesian product of a parameter grid.

    param_grid: dict
        {param_name: list_of_values}

    returns: list of dicts
        Each dict is one combination of parameters
    """
    keys = list(param_grid.keys())
    values = [param_grid[k] for k in keys]
    out = []
    for combo in itertools.product(*values):
        out.append(dict(zip(keys, combo)))
    return out


def stratified_bootstrap_sample(df, label_col="actual_positive", random_state=None):
    """
    Sample a dataframe with replacement, preserving class counts.

    Returns a dataframe with the same number of rows as the input, stratified on
    label_col.
    """
    rng = np.random.default_rng(random_state)
    sampled_parts = []

    for _, group in df.groupby(label_col, sort=False):
        sampled_idx = rng.choice(group.index.to_numpy(), size=len(group), replace=True)
        sampled_parts.append(df.loc[sampled_idx])

    boot_df = pd.concat(sampled_parts, axis=0)
    boot_df = boot_df.sample(frac=1, random_state=random_state).reset_index(drop=True)
    return boot_df


def resolve_cs_file(config):
    """
    Resolve CS file path from config.
    """
    if config.cs_filename:
        return Path(config.cs_out_dir) / config.cs_filename
    return Path(config.cs_out_dir) / f"{config.cs_run_id}.tsv"


# -----------------------------------------------------------------------------
# data preparation
# -----------------------------------------------------------------------------

def load_merged_for_config(config):
    """
    Load + collapse + merge CS and IC50 for a given config.
    """
    cs_file = resolve_cs_file(config)

    dr = load_drug_rankings(config.cs_out_dir, filename=cs_file.name)
    dr = add_pert_iname_to_cs(dr, config.lincs_metadata_path)

    dr = collapse_profiles_to_drug(
        dr,
        score_col=config.cs_score_col,
        drug_col=config.cs_drug_colname,
        how=config.cs_drug_collapse_method,
    )

    ic50, _ = load_IC50(
        config.ic50_file,
        config.disease,
        cell_line=config.cell_line_IC50,
        IC50_ONLY=config.ic50_only,
    )

    ic50 = collapse_profiles_to_drug(
        ic50,
        score_col=config.ic50_score_col,
        drug_col=config.ic50_drug_colname,
        how=config.ic50_drug_collapse_method,
    )

    merged = pd.merge(
        dr,
        ic50,
        left_on=config.cs_drug_colname,
        right_on=config.ic50_drug_colname,
        how="inner",
    )

    merged["actual_positive"] = (
        pd.to_numeric(merged[config.ic50_score_col], errors="coerce")
        <= config.ic50_threshold
    ).astype(int)

    merged = merged.dropna(subset=[config.cs_score_col, config.ic50_score_col])
    merged["drug_name"] = merged[config.cs_drug_colname].astype(str)

    return merged.reset_index(drop=True)


# -----------------------------------------------------------------------------
# bootstrap evaluation
# -----------------------------------------------------------------------------

def compute_metrics(eval_df, config):
    """
    Compute classification metrics on a dataframe.
    """
    classified_df, pr = classify_ic50_vs_cs(
        eval_df,
        score_col=config.cs_score_col,
        ic50_col=config.ic50_score_col,
        cs_threshold=config.cs_threshold,
        ic50_threshold=config.ic50_threshold,
    )

    y_true = classified_df["actual_positive"].astype(int)
    y_pred = classified_df["predicted_positive"].astype(int)
    f1 = f1_score(y_true, y_pred, zero_division=0)

    metrics = {
        "precision": pr["precision"],
        "recall": pr["recall"],
        "f1": f1,
        "tp": pr["tp"],
        "fp": pr["fp"],
        "fn": pr["fn"],
        "tn": pr["tn"],
    }
    return metrics, classified_df


def run_bootstrap_for_config(config, n_bootstraps=200, random_state=42):
    """
    Run stratified bootstrap evaluation for one config.

    returns
    -------
    full_classified : pd.DataFrame
        Full-dataset classified rows.
    boot_rows : pd.DataFrame
        One row per bootstrap replicate with precision / recall / F1.
    summary : dict
        Run-level summary metrics including full-dataset metrics and bootstrap
        mean / sd.
    """
    merged = load_merged_for_config(config)

    full_metrics, full_classified = compute_metrics(merged, config)
    full_metrics["bootstrap"] = "full_dataset"
    full_metrics["n_test"] = len(merged)

    boot_rows = []
    for i in range(n_bootstraps):
        boot_df = stratified_bootstrap_sample(merged,label_col="actual_positive",random_state=random_state + i,)
        metrics, _ = compute_metrics(boot_df, config)
        metrics["bootstrap"] = i
        metrics["n_test"] = len(boot_df)
        boot_rows.append(metrics)

    boot_rows = pd.DataFrame(boot_rows)

    summary = {
        "cv_run_name": TIMESTAMP,
        "config_name": config.name,
        "disease_run_id": disease_run_name,
        "drug_run_id": cell_line_run_name,
        "cs_run_id": config.cs_run_id,
        "n_rows_merged": len(merged),
        "n_unique_drugs": merged["drug_name"].nunique(),
        "n_positive": int(merged["actual_positive"].sum()),
        "n_negative": int((1 - merged["actual_positive"]).sum()),
        "n_bootstraps": n_bootstraps,
        "precision_full": full_metrics["precision"],
        "recall_full": full_metrics["recall"],
        "f1_full": full_metrics["f1"],
        "tp_full": full_metrics["tp"],
        "fp_full": full_metrics["fp"],
        "fn_full": full_metrics["fn"],
        "tn_full": full_metrics["tn"],
        "precision_mean": boot_rows["precision"].mean(),
        "precision_sd": boot_rows["precision"].std(ddof=1),
        "recall_mean": boot_rows["recall"].mean(),
        "recall_sd": boot_rows["recall"].std(ddof=1),
        "f1_mean": boot_rows["f1"].mean(),
        "f1_sd": boot_rows["f1"].std(ddof=1),
        "cs_threshold": config.cs_threshold,
        "ic50_threshold": config.ic50_threshold,
        "cs_drug_collapse_method": config.cs_drug_collapse_method,
        "ic50_drug_collapse_method": config.ic50_drug_collapse_method,
        "cell_line_IC50": config.cell_line_IC50,
        "ic50_only": config.ic50_only,
    }

    return full_classified, boot_rows, summary, full_metrics


# -----------------------------------------------------------------------------
# plotting
# -----------------------------------------------------------------------------

def plot_scatter_panel(classified_df, config, ax, metrics=None):
    """
    Plot one compact Fig.3-style scatter panel from full classified rows.

    metrics : dict or None
        Optional metrics dict. If provided, uses its precision/recall/f1 values
        instead of recomputing them inside the plotting function.
    """
    if classified_df is None or classified_df.empty:
        ax.set_axis_off()
        return

    x = pd.to_numeric(classified_df[config.cs_score_col], errors="coerce")
    y_raw = pd.to_numeric(classified_df[config.ic50_score_col], errors="coerce")
    ok = x.notna() & y_raw.notna() & (y_raw > 0)
    x = x[ok]
    y_raw = y_raw[ok]
    y = np.log10(y_raw)

    actual_pos = y_raw <= config.ic50_threshold
    edgecolors = np.where(actual_pos, "cornflowerblue", "lightcoral")

    ax.scatter(x, y, s=12, facecolors="none", edgecolors=edgecolors, linewidths=0.7)
    ax.axvline(config.cs_threshold, linestyle="--", linewidth=0.9, alpha=0.6)
    ax.axhline(np.log10(config.ic50_threshold), linestyle="--", linewidth=0.9, alpha=0.6)
    
    # label true positives
    tp_mask = actual_pos & (x <= config.cs_threshold)
    for xi, yi, name in zip(x[tp_mask], y[tp_mask], classified_df.loc[tp_mask, "drug_name"]):
        ax.text(xi,yi,name,fontsize=6,alpha=0.8,ha="left",va="bottom")

    if metrics is None:
        pred_pos = x <= config.cs_threshold
        tp = int((actual_pos & pred_pos).sum())
        fp = int((~actual_pos & pred_pos).sum())
        fn = int((actual_pos & ~pred_pos).sum())
        precision = tp / (tp + fp) if (tp + fp) else np.nan
        recall = tp / (tp + fn) if (tp + fn) else np.nan
        f1 = (
            2 * precision * recall / (precision + recall)
            if pd.notna(precision) and pd.notna(recall) and (precision + recall)
            else np.nan
        )
    else:
        precision = metrics.get("precision_full", metrics.get("precision", np.nan))
        recall = metrics.get("recall_full", metrics.get("recall", np.nan))
        f1 = metrics.get("f1_full", metrics.get("f1", np.nan))

    ax.set_title(config.name, fontsize=8)
    ax.text(
        0.02, 0.98,
        f"P={precision:.2f}  R={recall:.2f}  F1={f1:.2f}\n"
        f"N={len(x)}  cell={config.cell_line_IC50 if config.cell_line_IC50 is not None else 'ALL'}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
    )
    ax.set_xlabel("CS", fontsize=8)
    ax.set_ylabel("log10(IC50)", fontsize=8)
    ax.tick_params(labelsize=7)


def plot_config_matrix(all_full_classified, all_summaries, out_file, configs_by_name, metrics_by_name=None):
    """
    Plot a compact matrix of scatter panels, one per configuration.

    all_full_classified : dict
        {config_name: full_classified_df}
    all_summaries : pd.DataFrame or list of dict
        Summary table, ideally already sorted by preferred metric.
    out_file : str or Path
        Output image path.
    configs_by_name : dict
        {config_name: EvalConfig}
    metrics_by_name : dict or None
        Optional {config_name: metrics_dict} used to display already-computed
        full-dataset metrics in the plot text.
    """
    if isinstance(all_summaries, list):
        all_summaries = pd.DataFrame(all_summaries)

    all_summaries = all_summaries.sort_values("f1_mean", ascending=False).reset_index(drop=True)
    config_names = all_summaries["config_name"].tolist()

    n = len(config_names)
    if n == 0:
        return

    ncols = min(4, max(1, n))
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows), squeeze=False)
    axes = axes.ravel()

    for ax in axes:
        ax.set_axis_off()

    for ax, config_name in zip(axes, config_names):
        ax.set_axis_on()
        cfg = configs_by_name[config_name]
        classified_df = all_full_classified[config_name]
        metrics = None if metrics_by_name is None else metrics_by_name.get(config_name)

        plot_scatter_panel(classified_df, cfg, ax, metrics=metrics)

    fig.tight_layout()
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------

def make_default_config_grid():
    """
    Define hyperparameter grid.
    """
    run_id = resolve_cs_run_id(
        cs_runs_tsv=cs_log_filename,
        disease_run_id=disease_run_name,
        drug_run_id=cell_line_run_name,
        cs_on_LM=cs_on_LM,
        mith=cs_mith,
        selected_cs_run_id=selected_cs_run_id,
    )

    grid = {
        "cs_drug_collapse_method": ["best", "median"],
        "ic50_drug_collapse_method": ["best", "median"],
        "cs_threshold": [-1.5],
        "ic50_only": [True, False],
        "cell_line_IC50": [cell_line, None],
    }

    configs = []
    for i, params in enumerate(expand_grid(grid)):
        name = f"cfg_{i}"
        configs.append(
            EvalConfig(
                name=name,
                disease=DISEASE,
                cell_line=cell_line,
                cell_line_IC50=params["cell_line_IC50"],
                cs_out_dir=Path(CS_OUT),
                cs_run_id=run_id,
                cs_filename=f"{run_id}.tsv",
                cs_drug_collapse_method=params["cs_drug_collapse_method"],
                ic50_drug_collapse_method=params["ic50_drug_collapse_method"],
                cs_threshold=params["cs_threshold"],
                ic50_only=params["ic50_only"],
            )
        )

    return configs


#%%
if __name__ == "__main__":
    configs = make_default_config_grid()

    TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(LOGS_DIR) / "chembl_bootstrap_grid_search"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_summaries = []
    all_oof = {}
    configs_by_name = {cfg.name: cfg for cfg in configs}
    n_bootstraps = 200
    
    all_full_classified = {}
    all_full_metrics = {}
    configs_by_name = {cfg.name: cfg for cfg in configs}

    for cfg in configs:
        print("Running:", cfg.name)
        full_classified, boot_rows, summary, full_metrics = run_bootstrap_for_config(cfg, n_bootstraps=n_bootstraps)

        all_summaries.append(summary)
        all_oof[cfg.name] = full_classified
        all_full_classified[cfg.name] = full_classified
        all_full_metrics[cfg.name] = full_metrics
        print(pd.Series(summary))

        boot_rows.to_csv(out_dir / f"{DISEASE}_{TIMESTAMP}_{cfg.name}_bootstrap_metrics.tsv", sep="\t", index=False)
        full_classified.to_csv(out_dir / f"{DISEASE}_{TIMESTAMP}_{cfg.name}_full_dataset_classified.tsv", sep="\t", index=False)

    summary_df = pd.DataFrame(all_summaries).sort_values("f1_mean", ascending=False)
    summary_file = out_dir / f"{DISEASE}_{TIMESTAMP}_chembl_bootstrap_grid_search.tsv"
    summary_df.to_csv(summary_file, sep="\t", index=False)

    print("Saved bootstrap outputs to:", out_dir)
    matrix_file = IMG_DIR / f"{DISEASE}_{TIMESTAMP}_chembl_bootstrap_grid_search_matrix.png"
    plot_config_matrix(all_full_classified,summary_df,matrix_file, configs_by_name=configs_by_name, metrics_by_name=all_full_metrics,)
    # plot_config_matrix(all_oof, summary_df, matrix_file, configs_by_name=configs_by_name)
    print("Saved bootstrap images to:", IMG_DIR)

