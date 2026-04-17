# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 10:40:00 2026

Cross-validated hyperparameter search for ChEMBL IC50 vs connectivity score validation.

This script is intentionally built on top of final_step_correlations_IC50_CS.py.
It reuses its data-loading, collapsing, classification, and plotting helpers, then:

1. evaluates MANY fixed hyperparameter combinations under out-of-fold test evaluation
2. stores pooled out-of-fold predictions per configuration
3. summarizes precision / recall / F1 / Spearman across configurations
4. plots a matrix of fig3-like scatter panels, one panel per configuration
5. optionally performs nested CV to pick a single best configuration

Important scope note
--------------------
This script only tunes hyperparameters that can be evaluated from already-available
CS / IC50 tables. Hyperparameters that require recomputing the upstream CS itself
(e.g. LM vs all genes, DEG vs MITH, pathway-level vs gene-level signatures, different
CS runs) can still be included, but ONLY if the corresponding CS files / run IDs already
exist and are passed through the grid.

Recommended use
---------------
Start with the evaluation-only mode (outer CV, no nested selection). This gives an
honest out-of-fold estimate for each fixed configuration.

If later you want to *select one best configuration*, enable nested CV.

@Author: L-F-S
"""

from __future__ import annotations

import itertools
import math
import os
import sys
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, f1_score, precision_score, recall_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2] if len(HERE.parents) >= 3 else HERE.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from conf import CS_DIR, DATA_DIR, DISEASE,IC50_DRUG_COLLAPSE_METHOD,IC50_EFF_TH,\
    IC50_ONLY, IMG_DIR, LINCS_METADATA_PATH,LOGS_DIR, cell_line,cell_line_run_name,\
    chembl_val_log_filename,cs_filename,cs_log_filename,cs_mith,cs_on_LM,disease_run_name,\
    ic50_file, selected_cs_run_id
from logger import append_run_metadata 
from loader import load_IC50, load_drug_rankings
from final_step_correlations_IC50_CS import add_pert_iname_to_cs, classify_ic50_vs_cs,\
    collapse_profiles_to_drug, plot_binchen_fig3_style, resolve_cs_run_id

# -----------------------------------------------------------------------------
# config objects
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class EvalConfig:
    """Single hyperparameter configuration to evaluate."""

    # identity
    name: str

    # data source for CS
    disease: str
    cell_line: Optional[str]
    cs_out_dir: Path
    cs_run_id: Optional[str]
    cs_filename: Optional[str]

    # merge / collapse behavior
    cs_drug_colname: str = "pert_iname"
    ic50_drug_colname: str = "pert_iname"
    cs_score_col: str = "connectivity_score"
    ic50_score_col: str = "standard_value"
    cs_drug_collapse_method: str = "best"
    ic50_drug_collapse_method: str = "median"
    ic50_only: bool = True

    # decision thresholds
    cs_threshold: float = -1.5
    ic50_threshold: float = 10.0

    # labels / metadata
    lincs_metadata_path: Path = LINCS_METADATA_PATH
    ic50_file: Path = ic50_file

    # optional extra tags for reporting
    tags: Tuple[Tuple[str, Any], ...] = tuple()

    def to_dict(self) -> Dict[str, Any]:
        out = dict(self.__dict__)
        out["cs_out_dir"] = str(self.cs_out_dir)
        out["lincs_metadata_path"] = str(self.lincs_metadata_path)
        out["ic50_file"] = str(self.ic50_file)
        out["tags"] = dict(self.tags)
        return out


# -----------------------------------------------------------------------------
# utilities
# -----------------------------------------------------------------------------


def expand_grid(param_grid: Dict[str, Sequence[Any]]) -> List[Dict[str, Any]]:
    """Cartesian product of grid parameters."""
    keys = list(param_grid.keys())
    values = [param_grid[k] for k in keys]
    out = []
    for combo in itertools.product(*values):
        out.append(dict(zip(keys, combo)))
    return out



def make_default_config_grid() -> List[EvalConfig]:
    """
    Example grid.

    Adjust this to whatever you actually want to compare.
    For now it explores exactly the hyperparameters that are explicitly present in
    your validation slides / notes and are already supported by the current script:
      - CS drug collapse method
      - IC50 drug collapse method
      - CS threshold
      - IC50_ONLY
      - cell_line filter

    Any setting that changes the upstream CS file itself should be expressed by using
    different cs_run_id / cs_out_dir combinations.
    """
    resolved_run_id = resolve_cs_run_id(
        cs_runs_tsv=cs_log_filename,
        disease_run_id=disease_run_name,
        drug_run_id=cell_line_run_name,
        cs_on_LM=cs_on_LM,
        mith=cs_mith,
        selected_cs_run_id=selected_cs_run_id,
    )

    grid = {
        "cs_drug_collapse_method": ["best", "median", "mean"],
        "ic50_drug_collapse_method": ["median", "best"],
        "cs_threshold": [-1.95, -1.8, -1.65, -1.5, -1.35],
        "ic50_only": [True, False],
        "cell_line": [cell_line, None],
    }

    configs: List[EvalConfig] = []
    for i, params in enumerate(expand_grid(grid), start=1):
        name = (
            f"cfg_{i:03d}__cs-{params['cs_drug_collapse_method']}"
            f"__ic50-{params['ic50_drug_collapse_method']}"
            f"__thr-{params['cs_threshold']}"
            f"__IC50ONLY-{int(params['ic50_only'])}"
            f"__cell-{params['cell_line'] if params['cell_line'] is not None else 'ALL'}"
        )

        configs.append(
            EvalConfig(
                name=name,
                disease=DISEASE,
                cell_line=params["cell_line"],
                cs_out_dir=Path(CS_DIR) / "output" / DISEASE,
                cs_run_id=resolved_run_id,
                cs_filename=f"{resolved_run_id}.tsv",
                cs_drug_collapse_method=params["cs_drug_collapse_method"],
                ic50_drug_collapse_method=params["ic50_drug_collapse_method"],
                ic50_only=params["ic50_only"],
                cs_threshold=float(params["cs_threshold"]),
                ic50_threshold=float(IC50_EFF_TH),
                tags=tuple(sorted(params.items())),
            )
        )
    return configs



def pick_n_splits(y: Sequence[int], max_splits: int = 5, min_splits: int = 2) -> int:
    """
    Safe number of CV splits for small datasets.

    Stratified K-fold requires each class to have at least n_splits members.
    """
    y = np.asarray(y)
    _, counts = np.unique(y, return_counts=True)
    if len(counts) < 2:
        raise ValueError("Need both positive and negative classes for stratified CV.")
    n_splits = int(min(max_splits, counts.min()))
    if n_splits < min_splits:
        raise ValueError(
            f"Too few samples per class for CV. Class counts={counts.tolist()}, "
            f"need at least {min_splits}."
        )
    return n_splits



def ensure_numeric(df: pd.DataFrame, cols: Sequence[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out



def resolve_cs_file(config: EvalConfig) -> Path:
    if config.cs_filename is not None:
        return Path(config.cs_out_dir) / config.cs_filename
    if config.cs_run_id is not None:
        return Path(config.cs_out_dir) / f"{config.cs_run_id}.tsv"
    raise ValueError(f"Config '{config.name}' must define either cs_filename or cs_run_id.")


# -----------------------------------------------------------------------------
# data preparation per config
# -----------------------------------------------------------------------------


def load_merged_for_config(config: EvalConfig) -> pd.DataFrame:
    """Load, collapse, and merge CS and IC50 tables for a specific configuration."""
    cs_file = resolve_cs_file(config)

    dr_uncollapsed = load_drug_rankings(
        config.cs_out_dir,
        filename=cs_file.name,
    )
    dr_uncollapsed = add_pert_iname_to_cs(dr_uncollapsed, config.lincs_metadata_path)

    dr = collapse_profiles_to_drug(
        in_df=dr_uncollapsed,
        score_col=config.cs_score_col,
        drug_col=config.cs_drug_colname,
        how=config.cs_drug_collapse_method,
    )

    ic50_uncollapsed, _ = load_IC50(
        config.ic50_file,
        config.disease,
        cell_line=config.cell_line,
        IC50_ONLY=config.ic50_only,
    )

    ic50 = collapse_profiles_to_drug(
        in_df=ic50_uncollapsed,
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
        suffixes=("_dr", "_ic50"),
    )

    merged = ensure_numeric(merged, [config.cs_score_col, config.ic50_score_col])
    merged = merged.dropna(subset=[config.cs_score_col, config.ic50_score_col]).copy()
    merged["actual_positive"] = (merged[config.ic50_score_col] <= config.ic50_threshold).astype(int)
    merged["config_name"] = config.name
    merged["drug_name"] = merged[config.cs_drug_colname].astype(str)
    return merged.reset_index(drop=True)


# -----------------------------------------------------------------------------
# CV evaluation
# -----------------------------------------------------------------------------


def compute_fold_metrics(test_df: pd.DataFrame, config: EvalConfig) -> Dict[str, Any]:
    """Metrics on a single test fold using the configuration threshold."""
    classified_df, pr = classify_ic50_vs_cs(
        test_df,
        score_col=config.cs_score_col,
        ic50_col=config.ic50_score_col,
        cs_threshold=config.cs_threshold,
        ic50_threshold=config.ic50_threshold,
    )

    y_true = classified_df["actual_positive"].astype(int).to_numpy()
    y_pred = classified_df["predicted_positive"].astype(int).to_numpy()
    score = -classified_df[config.cs_score_col].astype(float).to_numpy()  # more positive = more likely active

    rho, rho_p = spearmanr(
        classified_df[config.cs_score_col].astype(float),
        np.log10(classified_df[config.ic50_score_col].astype(float)),
    )

    return {
        "n_test": int(len(classified_df)),
        "tp": int(pr["tp"]),
        "fp": int(pr["fp"]),
        "fn": int(pr["fn"]),
        "tn": int(pr["tn"]),
        "precision": float(pr["precision"]) if pd.notna(pr["precision"]) else np.nan,
        "recall": float(pr["recall"]) if pd.notna(pr["recall"]) else np.nan,
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "average_precision": float(average_precision_score(y_true, score)),
        "spearman_rho_log_ic50": float(rho) if pd.notna(rho) else np.nan,
        "spearman_p_log_ic50": float(rho_p) if pd.notna(rho_p) else np.nan,
    }, classified_df



def run_cv_for_config(
    config: EvalConfig,
    n_splits: int,
    n_repeats: int = 20,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    """
    Evaluate one fixed config by repeated stratified CV.

    Returns
    -------
    oof_df:
        pooled out-of-fold predictions / rows
    fold_df:
        per-fold metrics
    summary:
        summary metrics computed on pooled OOF data and fold means
    """
    merged = load_merged_for_config(config)
    y = merged["actual_positive"].astype(int).to_numpy()

    splitter = RepeatedStratifiedKFold(
        n_splits=n_splits,
        n_repeats=n_repeats,
        random_state=random_state,
    )

    oof_rows = []
    fold_rows = []

    for fold_id, (_, test_idx) in enumerate(splitter.split(merged, y), start=1):
        test_df = merged.iloc[test_idx].copy()
        fold_metrics, classified_df = compute_fold_metrics(test_df, config)
        fold_metrics["fold_id"] = fold_id
        fold_metrics["config_name"] = config.name
        fold_rows.append(fold_metrics)

        classified_df = classified_df.copy()
        classified_df["fold_id"] = fold_id
        classified_df["config_name"] = config.name
        oof_rows.append(classified_df)

    oof_df = pd.concat(oof_rows, axis=0, ignore_index=True)
    fold_df = pd.DataFrame(fold_rows)

    pooled_metrics_df, pooled_pr = classify_ic50_vs_cs(
        oof_df,
        score_col=config.cs_score_col,
        ic50_col=config.ic50_score_col,
        cs_threshold=config.cs_threshold,
        ic50_threshold=config.ic50_threshold,
    )

    y_true = pooled_metrics_df["actual_positive"].astype(int).to_numpy()
    y_pred = pooled_metrics_df["predicted_positive"].astype(int).to_numpy()
    pooled_score = -pooled_metrics_df[config.cs_score_col].astype(float).to_numpy()
    pooled_rho, pooled_rho_p = spearmanr(
        pooled_metrics_df[config.cs_score_col].astype(float),
        np.log10(pooled_metrics_df[config.ic50_score_col].astype(float)),
    )

    summary = {
        "config_name": config.name,
        "n_unique_drugs": int(merged[config.cs_drug_colname].nunique()),
        "n_rows_merged": int(len(merged)),
        "n_pos": int((merged["actual_positive"] == 1).sum()),
        "n_neg": int((merged["actual_positive"] == 0).sum()),
        "cs_drug_collapse_method": config.cs_drug_collapse_method,
        "ic50_drug_collapse_method": config.ic50_drug_collapse_method,
        "cs_threshold": float(config.cs_threshold),
        "ic50_threshold": float(config.ic50_threshold),
        "cell_line": config.cell_line if config.cell_line is not None else "ALL",
        "ic50_only": bool(config.ic50_only),
        "precision_oof": float(pooled_pr["precision"]) if pd.notna(pooled_pr["precision"]) else np.nan,
        "recall_oof": float(pooled_pr["recall"]) if pd.notna(pooled_pr["recall"]) else np.nan,
        "f1_oof": float(f1_score(y_true, y_pred, zero_division=0)),
        "average_precision_oof": float(average_precision_score(y_true, pooled_score)),
        "spearman_rho_log_ic50_oof": float(pooled_rho) if pd.notna(pooled_rho) else np.nan,
        "spearman_p_log_ic50_oof": float(pooled_rho_p) if pd.notna(pooled_rho_p) else np.nan,
        "precision_fold_mean": float(fold_df["precision"].mean()),
        "recall_fold_mean": float(fold_df["recall"].mean()),
        "f1_fold_mean": float(fold_df["f1"].mean()),
        "average_precision_fold_mean": float(fold_df["average_precision"].mean()),
        "precision_fold_sd": float(fold_df["precision"].std(ddof=1)),
        "recall_fold_sd": float(fold_df["recall"].std(ddof=1)),
        "f1_fold_sd": float(fold_df["f1"].std(ddof=1)),
        "average_precision_fold_sd": float(fold_df["average_precision"].std(ddof=1)),
    }
    return oof_df, fold_df, summary


# -----------------------------------------------------------------------------
# nested CV (optional)
# -----------------------------------------------------------------------------


def choose_best_config_on_train(
    train_indices: np.ndarray,
    base_data_by_config: Dict[str, pd.DataFrame],
    configs_by_name: Dict[str, EvalConfig],
    metric: str = "f1_oof",
    inner_random_state: int = 123,
) -> Tuple[str, pd.DataFrame]:
    """
    Select the best config using inner CV, *restricted to the outer train drugs*.

    This is only meaningful if all configs share the same set of drug names after merge.
    In practice they may differ because of IC50_ONLY / cell_line / collapse choices.
    So we reindex each config's table by drug_name intersection with the outer-train set.
    """
    inner_summaries = []

    # derive drug names from the first config available in base data
    first_name = next(iter(base_data_by_config))
    all_outer_train_drugs = set(base_data_by_config[first_name].iloc[train_indices]["drug_name"].tolist())

    for config_name, full_df in base_data_by_config.items():
        cfg = configs_by_name[config_name]
        sub_df = full_df[full_df["drug_name"].isin(all_outer_train_drugs)].copy().reset_index(drop=True)
        if sub_df["actual_positive"].nunique() < 2:
            continue

        try:
            inner_splits = pick_n_splits(sub_df["actual_positive"].to_numpy(), max_splits=4, min_splits=2)
        except ValueError:
            continue

        # temporary evaluation on train subset only
        y = sub_df["actual_positive"].astype(int).to_numpy()
        splitter = RepeatedStratifiedKFold(
            n_splits=inner_splits,
            n_repeats=5,
            random_state=inner_random_state,
        )
        fold_rows = []
        for _, test_idx in splitter.split(sub_df, y):
            test_df = sub_df.iloc[test_idx].copy()
            fold_metrics, _ = compute_fold_metrics(test_df, cfg)
            fold_rows.append(fold_metrics)

        fold_df = pd.DataFrame(fold_rows)
        inner_summary = {
            "config_name": config_name,
            "metric_value": float(fold_df[metric.replace("_oof", "").replace("_fold_mean", "")].mean())
            if metric.endswith("_fold_mean") else float(fold_df["f1"].mean()),
        }
        inner_summaries.append(inner_summary)

    inner_summary_df = pd.DataFrame(inner_summaries)
    if inner_summary_df.empty:
        raise ValueError("Nested CV could not evaluate any config on the outer train split.")

    best_row = inner_summary_df.sort_values("metric_value", ascending=False).iloc[0]
    return str(best_row["config_name"]), inner_summary_df


# -----------------------------------------------------------------------------
# plotting: matrix of fig3-like panels
# -----------------------------------------------------------------------------


def _plot_scatter_panel(ax: plt.Axes, df: pd.DataFrame, config: EvalConfig, summary: Dict[str, Any]) -> None:
    """Compact version of the left panel of plot_binchen_fig3_style, for grid plotting."""
    if df.empty:
        ax.set_axis_off()
        return

    x = pd.to_numeric(df[config.cs_score_col], errors="coerce").to_numpy(dtype=float)
    y = np.log10(pd.to_numeric(df[config.ic50_score_col], errors="coerce").to_numpy(dtype=float))
    positive_mask = pd.to_numeric(df[config.ic50_score_col], errors="coerce") <= config.ic50_threshold

    ax.scatter(
        x,
        y,
        s=8,
        facecolors="none",
        edgecolors=np.where(positive_mask, "cornflowerblue", "lightcoral"),
        lw=0.5,
    )
    ax.axvline(config.cs_threshold, linestyle="--", linewidth=1, color="red", alpha=0.2)
    ax.axhline(np.log10(config.ic50_threshold), linestyle="--", linewidth=1, color="red", alpha=0.2)

    txt = (
        f"P={summary['precision_oof']:.2f}  R={summary['recall_oof']:.2f}\n"
        f"F1={summary['f1_oof']:.2f}  AP={summary['average_precision_oof']:.2f}\n"
        f"N={summary['n_unique_drugs']}"
    )
    ax.text(
        0.03,
        0.97,
        txt,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75),
    )

    title = (
        f"{config.cs_drug_collapse_method} / {config.ic50_drug_collapse_method}\n"
        f"thr={config.cs_threshold:g} | cell={config.cell_line if config.cell_line is not None else 'ALL'}\n"
        f"IC50_ONLY={int(config.ic50_only)}"
    )
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("CS")
    ax.set_ylabel("log10(IC50)")



def plot_config_matrix(
    oof_by_config: Dict[str, pd.DataFrame],
    configs_by_name: Dict[str, EvalConfig],
    summary_df: pd.DataFrame,
    output_file: Path,
    ncols: int = 4,
    figsize_per_panel: Tuple[float, float] = (4.0, 3.2),
) -> None:
    """Matrix of compact fig3-like panels, one per config."""
    summary_lookup = summary_df.set_index("config_name").to_dict(orient="index")
    config_names = summary_df.sort_values("f1_oof", ascending=False)["config_name"].tolist()
    n_panels = len(config_names)
    nrows = math.ceil(n_panels / ncols)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
    )

    for ax in axes.flat:
        ax.set_axis_off()

    for ax, config_name in zip(axes.flat, config_names):
        ax.set_axis_on()
        _plot_scatter_panel(
            ax=ax,
            df=oof_by_config[config_name],
            config=configs_by_name[config_name],
            summary=summary_lookup[config_name],
        )

    fig.suptitle(
        f"Cross-validated ChEMBL validation matrix — {DISEASE}",
        fontsize=14,
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# full run
# -----------------------------------------------------------------------------


def run_grid_search(
    configs: Sequence[EvalConfig],
    out_dir: Path,
    n_repeats: int = 20,
    random_state: int = 42,
    nested_cv: bool = False,
    nested_metric: str = "f1_oof",
) -> Dict[str, Any]:
    """Main execution entrypoint."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # determine CV split count from the first config that loads cleanly
    first_df = load_merged_for_config(configs[0])
    n_splits = pick_n_splits(first_df["actual_positive"].to_numpy(), max_splits=5, min_splits=2)

    summaries = []
    oof_by_config: Dict[str, pd.DataFrame] = {}
    fold_by_config: Dict[str, pd.DataFrame] = {}
    configs_by_name = {cfg.name: cfg for cfg in configs}

    for cfg in configs:
        print(f"[grid] evaluating {cfg.name}")
        oof_df, fold_df, summary = run_cv_for_config(
            config=cfg,
            n_splits=n_splits,
            n_repeats=n_repeats,
            random_state=random_state,
        )
        summaries.append(summary)
        oof_by_config[cfg.name] = oof_df
        fold_by_config[cfg.name] = fold_df

        oof_df.to_csv(out_dir / f"{cfg.name}__oof.tsv", sep="\t", index=False)
        fold_df.to_csv(out_dir / f"{cfg.name}__folds.tsv", sep="\t", index=False)

    summary_df = pd.DataFrame(summaries).sort_values("f1_oof", ascending=False).reset_index(drop=True)
    summary_file = out_dir / "grid_search_summary.tsv"
    summary_df.to_csv(summary_file, sep="\t", index=False)

    matrix_file = out_dir / f"{DISEASE}_cv_grid_matrix.png"
    plot_config_matrix(
        oof_by_config=oof_by_config,
        configs_by_name=configs_by_name,
        summary_df=summary_df,
        output_file=matrix_file,
    )

    nested_results = None
    if nested_cv:
        print("[nested] running nested CV selection")
        # nested CV uses only the first config's drugs to define outer splits.
        base_df = load_merged_for_config(configs[0])
        y = base_df["actual_positive"].astype(int).to_numpy()
        outer_splits = pick_n_splits(y, max_splits=5, min_splits=2)
        outer_splitter = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=random_state)

        base_data_by_config = {cfg.name: load_merged_for_config(cfg) for cfg in configs}
        chosen_rows = []
        for outer_fold, (train_idx, test_idx) in enumerate(outer_splitter.split(base_df, y), start=1):
            best_name, inner_df = choose_best_config_on_train(
                train_indices=train_idx,
                base_data_by_config=base_data_by_config,
                configs_by_name=configs_by_name,
                metric=nested_metric,
                inner_random_state=random_state + outer_fold,
            )

            chosen_cfg = configs_by_name[best_name]
            test_drugs = set(base_df.iloc[test_idx]["drug_name"].tolist())
            test_df = base_data_by_config[best_name]
            test_df = test_df[test_df["drug_name"].isin(test_drugs)].copy().reset_index(drop=True)
            fold_metrics, _ = compute_fold_metrics(test_df, chosen_cfg)
            fold_metrics["outer_fold"] = outer_fold
            fold_metrics["selected_config"] = best_name
            chosen_rows.append(fold_metrics)

            inner_df.to_csv(out_dir / f"nested_outerfold_{outer_fold:02d}__inner_summary.tsv", sep="\t", index=False)

        nested_results = pd.DataFrame(chosen_rows)
        nested_results.to_csv(out_dir / "nested_cv_outer_results.tsv", sep="\t", index=False)

    log_payload = {
        "datetime": datetime.now().isoformat(),
        "disease": DISEASE,
        "n_configs": len(configs),
        "n_splits": n_splits,
        "n_repeats": n_repeats,
        "nested_cv": bool(nested_cv),
        "summary_file": str(summary_file),
        "matrix_file": str(matrix_file),
        "best_config_by_f1_oof": summary_df.iloc[0]["config_name"],
        "best_f1_oof": float(summary_df.iloc[0]["f1_oof"]),
        "best_precision_oof": float(summary_df.iloc[0]["precision_oof"]),
        "best_recall_oof": float(summary_df.iloc[0]["recall_oof"]),
    }
    append_run_metadata(chembl_val_log_filename, log_payload)

    return {
        "summary_df": summary_df,
        "summary_file": summary_file,
        "matrix_file": matrix_file,
        "nested_results": nested_results,
        "oof_by_config": oof_by_config,
        "fold_by_config": fold_by_config,
    }


# -----------------------------------------------------------------------------
# command-line execution
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    results_dir = Path(LOGS_DIR) / "chembl_cv_grid_search" / datetime.now().strftime("%Y%m%d_%H%M%S")
    configs = make_default_config_grid()
    run_grid_search(
        configs=configs,
        out_dir=results_dir,
        n_repeats=20,
        random_state=42,
        nested_cv=False,   # switch to True only when you really want model selection
        nested_metric="f1_oof",
    )
    print(f"Saved results to: {results_dir}")