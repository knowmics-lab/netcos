# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 10:40:00 2026

Bootstrap hyperparameter search for ChEMBL IC50 vs connectivity score validation.

Mirrors the structure of call_cs_batch_for_different_confs.py:

  1. import the shared grid (upstream pre-CS knobs + downstream validation knobs)
     from validations.chembl.grid
  2. enumerate upstream x downstream combos
  3. drop combos that violate conf.validate_hyperparameters
  4. drop upstream combos that have no matching row in logs/cs_runs.tsv
     (i.e. the CS file does not exist locally)
  5. for each remaining combo, call run_correlation_for_conf(...) from
     final_step_correlations_IC50_CS.py to get the merged classified dataframe
     and full-dataset precision/recall metrics
  6. run a stratified bootstrap on the merged dataframe to get
     precision / recall / F1 mean and sd
  7. save per-combo bootstrap TSV, per-combo full-dataset classified TSV, and
     a summary TSV under logs/chembl_bootstrap_grid_search/

Important scope note
--------------------
This script does not retrain anything; it evaluates fixed configurations.
The matrix plot (plot_config_matrix) is intentionally NOT called here — the
plotting helpers are kept in this file but the user will plot from the saved
TSVs and CS files later.

@Author: L-F-S
"""

from __future__ import annotations

import itertools
import math
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

from conf import (
    CS_DIR, CS_OUT, DATA_DIR, IMG_DIR, LINCS_METADATA_PATH, LOGS_DIR,
    chembl_val_log_filename, cs_log_filename, diseases_of, ic50_file,
    validate_hyperparameters,
)
from validations.chembl.final_step_correlations_IC50_CS import (
    classify_ic50_vs_cs, run_correlation_for_conf,
)
from validations.chembl.grid import (
    cell_lines, landmark_diseases, landmark_drugs,
    CS_METHODs, cs_miths, cs_on_LMs, CS_ON_PATHWAYSs,
    cs_drug_collapse_methods, ic50_drug_collapse_methods,
    ic50_onlys, cell_line_IC50_options, cs_thresholds,
)


# -----------------------------------------------------------------------------
# config object (one row of the per-combo summary)
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class EvalConfig:
    """Container for one (upstream + downstream) hyperparameter combination."""
    name: str
    # pre-CS
    cell_line: str
    disease: str
    landmark_disease: bool
    landmark_drug: bool
    CS_METHOD: str
    cs_mith: int
    cs_on_LM: int
    CS_ON_PATHWAYS: bool
    # validation
    cs_drug_collapse_method: str
    ic50_drug_collapse_method: str
    ic50_only: bool
    cell_line_IC50: object  # str or None
    cs_threshold: float = -1.5
    ic50_threshold: float = 10.0


# -----------------------------------------------------------------------------
# grid utilities
# -----------------------------------------------------------------------------

def _to_bool(v):
    """Coerce bool/int/float/str values from cs_runs.tsv into a Python bool."""
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return bool(int(v))
    s = str(v).strip().lower()
    return s in ('true', '1', '1.0')


def has_cs_run(log_df, cell_line, landmark_disease, landmark_drug,
               CS_METHOD, cs_mith, cs_on_LM, CS_ON_PATHWAYS):
    """
    Inverse of call_cs_batch_for_different_confs.already_logged: True iff the
    CS row exists in cs_runs.tsv (i.e. the CS file is available to validate).
    """
    if log_df.empty:
        return False
    disease_run_id = diseases_of[cell_line] + ('_LM' if landmark_disease else '')
    drug_run_id    = cell_line              + ('_LM' if landmark_drug    else '')
    log_pw_bool = log_df['CS_ON_PATHWAYS'].apply(_to_bool)
    match = (
        (log_df['disease_run_id'] == disease_run_id) &
        (log_df['drug_run_id']    == drug_run_id) &
        (log_df['mith'].astype(int)        == int(cs_mith)) &
        (log_df['cs_on_LM'].astype(int)    == int(cs_on_LM)) &
        (log_pw_bool                       == bool(CS_ON_PATHWAYS)) &
        (log_df['CS_METHOD']               == CS_METHOD)
    )
    return bool(match.any())


def format_combo(combo):
    cell_line, ld_disease, ld_drug, method, mith, on_lm, on_pw = combo
    return (f"cell_line={cell_line} landmark_disease={ld_disease} "
            f"landmark_drug={ld_drug} CS_METHOD={method} "
            f"cs_mith={mith} cs_on_LM={on_lm} CS_ON_PATHWAYS={on_pw}")


# -----------------------------------------------------------------------------
# bootstrap utilities
# -----------------------------------------------------------------------------

def stratified_bootstrap_sample(df, label_col="actual_positive", random_state=None):
    """
    Sample a dataframe with replacement, preserving class counts.
    Returns a dataframe with the same number of rows as the input, stratified
    on label_col.
    """
    rng = np.random.default_rng(random_state)
    sampled_parts = []
    for _, group in df.groupby(label_col, sort=False):
        sampled_idx = rng.choice(group.index.to_numpy(), size=len(group), replace=True)
        sampled_parts.append(df.loc[sampled_idx])
    boot_df = pd.concat(sampled_parts, axis=0)
    boot_df = boot_df.sample(frac=1, random_state=random_state).reset_index(drop=True)
    return boot_df


def compute_metrics_on_df(eval_df, cs_score_col, ic50_score_col, cs_threshold, ic50_threshold):
    """Compute classification metrics on a dataframe (one bootstrap replicate)."""
    classified_df, pr = classify_ic50_vs_cs(
        eval_df,
        score_col=cs_score_col,
        ic50_col=ic50_score_col,
        cs_threshold=cs_threshold,
        ic50_threshold=ic50_threshold,
    )
    y_true = classified_df["actual_positive"].astype(int)
    y_pred = classified_df["predicted_positive"].astype(int)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    return {
        "precision": pr["precision"],
        "recall": pr["recall"],
        "f1": f1,
        "tp": pr["tp"],
        "fp": pr["fp"],
        "fn": pr["fn"],
        "tn": pr["tn"],
    }


def run_bootstrap(merged_classified, cfg, n_bootstraps=200, random_state=42):
    """
    Run stratified bootstrap on a pre-classified merged dataframe.

    Returns
    -------
    boot_rows : pd.DataFrame    one row per bootstrap replicate
    boot_summary : dict         mean / sd of precision / recall / F1
    """
    # classify_ic50_vs_cs adds 'actual_positive' / 'predicted_positive' columns;
    # we resample using actual_positive as the stratifying label
    boot_rows = []
    for i in range(n_bootstraps):
        boot_df = stratified_bootstrap_sample(
            merged_classified,
            label_col="actual_positive",
            random_state=random_state + i,
        )
        m = compute_metrics_on_df(
            boot_df,
            cs_score_col="connectivity_score",
            ic50_score_col="standard_value",
            cs_threshold=cfg.cs_threshold,
            ic50_threshold=cfg.ic50_threshold,
        )
        m["bootstrap"] = i
        m["n_test"] = len(boot_df)
        boot_rows.append(m)

    boot_rows = pd.DataFrame(boot_rows)
    boot_summary = {
        "n_bootstraps": n_bootstraps,
        "precision_mean": boot_rows["precision"].mean(),
        "precision_sd": boot_rows["precision"].std(ddof=1),
        "recall_mean": boot_rows["recall"].mean(),
        "recall_sd": boot_rows["recall"].std(ddof=1),
        "f1_mean": boot_rows["f1"].mean(),
        "f1_sd": boot_rows["f1"].std(ddof=1),
    }
    return boot_rows, boot_summary


# -----------------------------------------------------------------------------
# plotting helpers (kept for later use; not called from main)
# -----------------------------------------------------------------------------

def plot_scatter_panel(classified_df, config, ax, metrics=None):
    """Compact Fig.3-style scatter panel for the matrix plot."""
    if classified_df is None or classified_df.empty:
        ax.set_axis_off()
        return

    x = pd.to_numeric(classified_df["connectivity_score"], errors="coerce")
    y_raw = pd.to_numeric(classified_df["standard_value"], errors="coerce")
    ok = x.notna() & y_raw.notna() & (y_raw > 0)
    x = x[ok]
    y_raw = y_raw[ok]
    y = np.log10(y_raw)

    actual_pos = y_raw <= config.ic50_threshold
    edgecolors = np.where(actual_pos, "cornflowerblue", "lightcoral")

    ax.scatter(x, y, s=12, facecolors="none", edgecolors=edgecolors, linewidths=0.7)
    ax.axvline(config.cs_threshold, linestyle="--", linewidth=0.9, alpha=0.6)
    ax.axhline(np.log10(config.ic50_threshold), linestyle="--", linewidth=0.9, alpha=0.6)

    tp_mask = actual_pos & (x <= config.cs_threshold)
    drug_names = classified_df.loc[ok, "pert_iname"] if "pert_iname" in classified_df.columns else None
    if drug_names is not None:
        for xi, yi, name in zip(x[tp_mask], y[tp_mask], drug_names[tp_mask]):
            ax.text(xi, yi, name, fontsize=6, alpha=0.8, ha="left", va="bottom")

    if metrics is None:
        pred_pos = x <= config.cs_threshold
        tp = int((actual_pos & pred_pos).sum())
        fp = int((~actual_pos & pred_pos).sum())
        fn = int((actual_pos & ~pred_pos).sum())
        precision = tp / (tp + fp) if (tp + fp) else np.nan
        recall = tp / (tp + fn) if (tp + fn) else np.nan
        f1 = (2 * precision * recall / (precision + recall)
              if pd.notna(precision) and pd.notna(recall) and (precision + recall)
              else np.nan)
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
        va="top", ha="left", fontsize=7,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
    )
    ax.set_xlabel("CS", fontsize=8)
    ax.set_ylabel("log10(IC50)", fontsize=8)
    ax.tick_params(labelsize=7)


def plot_config_matrix(all_full_classified, all_summaries, out_file, configs_by_name, metrics_by_name=None):
    """
    Compact matrix of scatter panels, one per configuration.
    Not called from main — kept as a helper for offline plotting from the
    saved per-combo TSVs.
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

def main(n_bootstraps=200, ic50_threshold=10.0):
    TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(LOGS_DIR) / "chembl_bootstrap_grid_search"
    out_dir.mkdir(parents=True, exist_ok=True)

    log_df = (pd.read_csv(cs_log_filename, sep='\t')
              if cs_log_filename.exists() and cs_log_filename.stat().st_size > 0
              else pd.DataFrame())

    # ---- 1. enumerate the full hyperparameter grid in one product ----
    all_combos = list(itertools.product(
        cell_lines, landmark_diseases, landmark_drugs,
        CS_METHODs, cs_miths, cs_on_LMs, CS_ON_PATHWAYSs,
        cs_drug_collapse_methods, ic50_drug_collapse_methods,
        ic50_onlys, cell_line_IC50_options, cs_thresholds))

    # ---- 2. filter and build EvalConfigs in one pass ----
    configs = []
    n_invalid = 0
    no_cs_row_upstream = set()
    for combo in all_combos:
        (cell_line, ld_dis, ld_drug, method, mith, on_lm, on_pw,
         cs_collapse, ic50_collapse, ic50_only, cl_ic50_opt, cs_thr) = combo

        try:
            validate_hyperparameters(mith, on_lm, on_pw,
                                     CS_METHOD=method, cs_threshold=cs_thr)
        except ValueError:
            n_invalid += 1
            continue

        upstream = (cell_line, ld_dis, ld_drug, method, mith, on_lm, on_pw)
        if not has_cs_run(log_df, *upstream):
            no_cs_row_upstream.add(upstream)
            continue

        cl_ic50 = cell_line if cl_ic50_opt == 'per_cell_line' else cl_ic50_opt
        idx = len(configs)
        name = (
            f"cfg_{idx:04d}_{cell_line}"
            f"_{'LMd' if ld_dis else 'NLd'}"
            f"_{method.replace('_', '')}"
            f"_mith{mith}_lm{on_lm}_pw{int(on_pw)}"
            f"_cs{cs_collapse}_ic{ic50_collapse}"
            f"_only{int(ic50_only)}_cl{cl_ic50 or 'ALL'}"
            f"_th{cs_thr}"
        )
        configs.append(EvalConfig(
            name=name,
            cell_line=cell_line,
            disease=diseases_of[cell_line],
            landmark_disease=ld_dis,
            landmark_drug=ld_drug,
            CS_METHOD=method,
            cs_mith=mith,
            cs_on_LM=on_lm,
            CS_ON_PATHWAYS=on_pw,
            cs_drug_collapse_method=cs_collapse,
            ic50_drug_collapse_method=ic50_collapse,
            ic50_only=ic50_only,
            cell_line_IC50=cl_ic50,
            cs_threshold=cs_thr,
            ic50_threshold=ic50_threshold,
        ))

    n_no_cs_row = len(all_combos) - len(configs) - n_invalid

    print(f"Total combos                     : {len(all_combos)}")
    print(f"  invalid (hyperparameter check) : {n_invalid}")
    print(f"  missing CS row                 : {n_no_cs_row}")
    print(f"  evaluating                     : {len(configs)}")

    if no_cs_row_upstream:
        print("\nUpstream tuples with no CS row in cs_runs.tsv (skipped):")
        for u in sorted(no_cs_row_upstream):
            print(f"  - {format_combo(u)}")

    if not configs:
        print("\nNothing to evaluate. Exiting.")
        return

    # ---- 4. evaluate each config ----
    all_summaries = []                # one row per combo, for the grid-search summary TSV
    all_full_classified = {}
    all_full_metrics = {}
    chembl_log_rows = []              # batched at the end into chembl_val_log_filename
    configs_by_name = {cfg.name: cfg for cfg in configs}

    for i, cfg in enumerate(configs, 1):
        print(f"\n=== [{i}/{len(configs)}] {cfg.name} ===")
        try:
            classified_df, pr_metrics, summary = run_correlation_for_conf(
                disease=cfg.disease,
                cell_line=cfg.cell_line,
                landmark_disease=cfg.landmark_disease,
                landmark_drug=cfg.landmark_drug,
                CS_METHOD=cfg.CS_METHOD,
                cs_mith=cfg.cs_mith,
                cs_on_LM=cfg.cs_on_LM,
                CS_ON_PATHWAYS=cfg.CS_ON_PATHWAYS,
                cs_drug_collapse_method=cfg.cs_drug_collapse_method,
                ic50_drug_collapse_method=cfg.ic50_drug_collapse_method,
                ic50_only=cfg.ic50_only,
                cell_line_IC50=cfg.cell_line_IC50,
                cs_threshold=cfg.cs_threshold,
                ic50_threshold=cfg.ic50_threshold,
                compute_sd5=False,
                plot=False,
                log_to_tsv=False,           # batched at the end (see below)
                verbose=False,
            )
        except Exception as e:
            print(f"!!! Config failed ({cfg.name}): {e}", file=sys.stderr)
            continue

        # bootstrap on the merged classified dataframe
        boot_rows, boot_summary = run_bootstrap(
            classified_df, cfg, n_bootstraps=n_bootstraps,
        )
        # merge full-dataset metrics into the per-combo summary
        full_metrics = {
            "precision_full": pr_metrics["precision"],
            "recall_full": pr_metrics["recall"],
            "f1_full": (
                2 * pr_metrics["precision"] * pr_metrics["recall"]
                / (pr_metrics["precision"] + pr_metrics["recall"])
                if pd.notna(pr_metrics["precision"]) and pd.notna(pr_metrics["recall"])
                and (pr_metrics["precision"] + pr_metrics["recall"]) > 0
                else np.nan
            ),
            "tp_full": pr_metrics["tp"],
            "fp_full": pr_metrics["fp"],
            "fn_full": pr_metrics["fn"],
            "tn_full": pr_metrics["tn"],
        }

        combo_summary = {
            "cv_run_name": TIMESTAMP,
            "config_name": cfg.name,
            **{k: getattr(cfg, k) for k in (
                "cell_line", "landmark_disease", "landmark_drug",
                "CS_METHOD", "cs_mith", "cs_on_LM", "CS_ON_PATHWAYS",
                "cs_drug_collapse_method", "ic50_drug_collapse_method",
                "ic50_only", "cell_line_IC50", "cs_threshold", "ic50_threshold",
            )},
            "cs_run_id": summary["cs_run_id"],
            "n_rows_merged": summary["n_merged_rows"],
            "n_unique_drugs": summary["n_unique_drugs"],
            **full_metrics,
            **boot_summary,
        }
        all_summaries.append(combo_summary)
        all_full_classified[cfg.name] = classified_df
        all_full_metrics[cfg.name] = full_metrics
        chembl_log_rows.append(summary)

        # save per-combo bootstrap rows + full-dataset classified rows
        disease_id = diseases_of[cfg.cell_line] + ('_LM' if cfg.landmark_disease else '')
        boot_rows.to_csv(
            out_dir / f"{disease_id}_{TIMESTAMP}_{cfg.name}_bootstrap_metrics.tsv",
            sep="\t", index=False,
        )
        classified_df.to_csv(
            out_dir / f"{disease_id}_{TIMESTAMP}_{cfg.name}_full_dataset_classified.tsv",
            sep="\t", index=False,
        )

    # ---- 5. write the grid-search summary TSV ----
    summary_df = pd.DataFrame(all_summaries).sort_values("f1_mean", ascending=False)
    summary_file = out_dir / f"{TIMESTAMP}_chembl_bootstrap_grid_search.tsv"
    summary_df.to_csv(summary_file, sep="\t", index=False)
    print("\nSaved bootstrap outputs to:", out_dir)
    print("Summary TSV:", summary_file)

    # ---- 6. batch-append per-combo summaries to chembl_val_log_filename ----
    # One read + one write at the end avoids the O(N^2) read/rewrite pattern
    # that per-iteration logging would have caused.
    if chembl_log_rows:
        new_rows_df = pd.DataFrame(chembl_log_rows)
        if (chembl_val_log_filename.exists()
                and chembl_val_log_filename.stat().st_size > 0):
            old_df = pd.read_csv(chembl_val_log_filename, sep='\t')
            out_df = pd.concat([old_df, new_rows_df], ignore_index=True)
        else:
            out_df = new_rows_df
        chembl_val_log_filename.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(chembl_val_log_filename, sep='\t', index=False)
        print(f"Appended {len(new_rows_df)} rows to {chembl_val_log_filename}")

    # ---- 7. print top-10 by precision and by recall ----
    if not summary_df.empty:
        cols = [
            'config_name', 'cell_line', 'landmark_disease', 'CS_METHOD',
            'cs_mith', 'cs_on_LM', 'CS_ON_PATHWAYS',
            'cs_drug_collapse_method', 'ic50_drug_collapse_method',
            'cs_threshold', 'precision_mean', 'recall_mean', 'f1_mean',
            'precision_full', 'recall_full', 'n_rows_merged',
        ]
        cols = [c for c in cols if c in summary_df.columns]
        with pd.option_context('display.max_colwidth', None,
                               'display.width', 200,
                               'display.float_format', '{:.3f}'.format):
            print("\n=== Top 10 by precision_mean ===")
            print(summary_df.sort_values('precision_mean', ascending=False)
                            .head(10)[cols].to_string(index=False))
            print("\n=== Top 10 by recall_mean ===")
            print(summary_df.sort_values('recall_mean', ascending=False)
                            .head(10)[cols].to_string(index=False))

    print("\nSkipping matrix plot — will be generated later from the saved TSVs and CS files.")


if __name__ == "__main__":
    main()
