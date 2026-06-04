#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026-05-13

@author: L-F-S

sRGES replication check against Bin Chen 2017 reference data.

Three scopes are checked:

A) `--scope full`           [BinChen RGES -> my sRGES vs BinChen sRGES]
   Input  : data/BinChen2017/data/data/<DISEASE>/all_lincs_score.csv
   Target : data/BinChen2017/data/data/<DISEASE>/lincs_cancer_sRGES.csv
   Tests *only the sRGES aggregation algorithm*. Per-instance RGES values
   come from Bin Chen's already-computed CSV, so any disagreement isolates
   bugs in the aggregation step (build_reference_diff / getsRGES port /
   per-drug mean).

B) `--scope per_cell_line`  [BinChen RGES, restricted -> my sRGES vs BinChen sRGES]
   Input  : same all_lincs_score.csv, but filtered to one cell_id.
   Target : data/BinChen2017/data/data/<DISEASE>/rges_ic50_normalized.csv
            (per-drug sRGES + IC50, drugs with IC50 only).
   Figure 3 single-cell-line variant (LIHC->HEPG2, BRCA->MCF7, COAD->HT29).

C) `--scope full_netcos`    [NetCoS RGES -> my sRGES vs BinChen sRGES]
   Input  : a NetCoS CS .tsv (or several, concatenated across cell lines)
            from connectivity_score/output/<DISEASE>/, enriched with the
            LINCS metadata (pert_iname / cell_id / pert_dose / pert_time /
            is_gold) via lincs_sig_info_new.csv.
   Target : same lincs_cancer_sRGES.csv.
   This is the end-to-end test the project actually cares about: it stacks
   the per-signature RGES bit-exactness already proven by
   compare_rges_per_disease.py on top of the aggregation algorithm.

Everything that is NOT specific to this comparison driver lives elsewhere:
  - loaders            -> src/loader.py
  - default parameters -> src/conf.py
  - sRGES algorithm    -> src/sRGES.py

Usage
-----
    # Aggregation-only test, LIHC
    python -m src.validations.binchen2017.test_sRGES_replication \
        --disease LIHC --scope full --use-cor

    # End-to-end NetCoS test, LIHC, auto-picking the canonical CS run
    python -m src.validations.binchen2017.test_sRGES_replication \
        --disease LIHC --scope full_netcos --use-cor

    # End-to-end NetCoS test, explicit per-cell-line CS files
    python -m src.validations.binchen2017.test_sRGES_replication \
        --disease LIHC --scope full_netcos --use-cor \
        --netcos-cs-file connectivity_score/output/LIHC/<run_HEPG2>.tsv \
        --netcos-cs-file connectivity_score/output/LIHC/<run_MCF7>.tsv \
        --netcos-cs-file connectivity_score/output/LIHC/<run_HT29>.tsv

    # Per-cell-line BinChen test (Fig 3 cell line for the disease)
    python -m src.validations.binchen2017.test_sRGES_replication \
        --disease LIHC --scope per_cell_line --use-cor

    # All diseases x all three scopes
    python -m src.validations.binchen2017.test_sRGES_replication --all --use-cor
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Allow running as a script from anywhere within the repo:
THIS_DIR = Path(__file__).resolve().parent
SRC_DIR = THIS_DIR.parent.parent  # .../src
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import conf  # noqa: E402
from sRGES import compute_sRGES  # noqa: E402
from loader import (  # noqa: E402
    cs_out_dir_for,
    pick_canonical_cs_run,
    load_netcos_cs,
    load_netcos_cs_combined,
    load_binchen_all_lincs_score,
    load_reference_full_sRGES,
    load_reference_fig3_sRGES,
    cor_weights_for_disease,
)


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

def compare(my_df, ref_df, *, key_col="pert_iname",
            my_score_col="sRGES", ref_score_col="sRGES"):
    """Join my_df and ref_df on key_col, compute Pearson/Spearman + topK overlap."""
    merged = my_df.merge(ref_df, on=key_col, suffixes=("_mine", "_ref"))
    if len(merged) < 3:
        return merged, {
            "n_joined": len(merged),
            "pearson_r": np.nan, "pearson_p": np.nan,
            "spearman_rho": np.nan, "spearman_p": np.nan,
            "rmse": np.nan, "mae": np.nan,
        }
    mine = (merged[my_score_col + "_mine"]
            if my_score_col + "_mine" in merged.columns
            else merged[my_score_col])
    refv = (merged[ref_score_col + "_ref"]
            if ref_score_col + "_ref" in merged.columns
            else merged[ref_score_col])

    pear = stats.pearsonr(mine, refv)
    spe = stats.spearmanr(mine, refv)
    diff = (mine - refv).to_numpy()

    metrics = {
        "n_joined": len(merged),
        "pearson_r": float(pear.statistic if hasattr(pear, "statistic") else pear[0]),
        "pearson_p": float(pear.pvalue if hasattr(pear, "pvalue") else pear[1]),
        "spearman_rho": float(spe.statistic if hasattr(spe, "statistic") else spe.correlation),
        "spearman_p": float(spe.pvalue if hasattr(spe, "pvalue") else spe.pvalue),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "mae": float(np.mean(np.abs(diff))),
    }
    # top-K overlap (most-negative sRGES)
    for k in (10, 50, 100):
        if len(merged) >= k:
            top_mine = set(merged.nsmallest(k, mine.name).index) if hasattr(mine, "name") else set()
            top_ref = set(merged.nsmallest(k, refv.name).index) if hasattr(refv, "name") else set()
            metrics[f"top{k}_overlap"] = len(top_mine & top_ref) / k if k else np.nan
    return merged, metrics


def plot_scatter(merged, *, my_col, ref_col, title, out_path):
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    x = merged[ref_col].to_numpy()
    y = merged[my_col].to_numpy()
    ax.scatter(x, y, s=10, alpha=0.5, edgecolor="none")
    lo = float(min(np.nanmin(x), np.nanmin(y)))
    hi = float(max(np.nanmax(x), np.nanmax(y)))
    ax.plot([lo, hi], [lo, hi], "k--", linewidth=1, alpha=0.7, label="y=x")
    ax.set_xlabel(f"BinChen ref ({ref_col})")
    ax.set_ylabel(f"NetCoS sRGES ({my_col})")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def _print_metrics(header, metrics):
    print(f"\n=== {header} ===")
    print(f"  n joined drugs           : {metrics['n_joined']}")
    print(f"  Pearson  r vs BinChen    : {metrics['pearson_r']:.6f}  "
          f"(p={metrics['pearson_p']:.2e})")
    print(f"  Spearman rho vs BinChen  : {metrics['spearman_rho']:.6f}  "
          f"(p={metrics['spearman_p']:.2e})")
    print(f"  RMSE / MAE               : {metrics['rmse']:.6f} / "
          f"{metrics['mae']:.6f}")
    for k in (10, 50, 100):
        if f"top{k}_overlap" in metrics:
            print(f"  top-{k} overlap (low)     : {metrics[f'top{k}_overlap']:.3f}")


def _maybe_save(merged, metrics, *, out_dir, disease, tag, title):
    if out_dir is None:
        return
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_dir / f"{disease}__{tag}__joined.tsv",
                  sep="\t", index=False)
    pd.DataFrame([metrics]).to_csv(
        out_dir / f"{disease}__{tag}__summary.tsv", sep="\t", index=False
    )
    plot_scatter(
        merged, my_col="sRGES_mine", ref_col="sRGES_ref",
        title=title,
        out_path=out_dir / f"{disease}__{tag}__scatter.png",
    )


# ---------------------------------------------------------------------------
# Main per-disease drivers
# ---------------------------------------------------------------------------

def run_full_sRGES(disease, *, use_cor=False,
                   srges_diff_mode=conf.srges_diff_mode,
                   out_dir=None):
    """Scope A: BinChen RGES -> my sRGES vs BinChen sRGES (aggregation only)."""
    inp = load_binchen_all_lincs_score(disease)
    ref = load_reference_full_sRGES(disease)

    weights = cor_weights_for_disease(disease) if use_cor else None
    my = compute_sRGES(
        inp, score_col="RGES",
        cell_lines=None,                # full -> all cell lines
        cell_line_weights=weights,
        apply_cor_weight=use_cor,
        filter_is_gold=False,
        srges_diff_mode=srges_diff_mode,
    )

    merged, metrics = compare(my, ref, key_col="pert_iname")

    tag = f"full__corweighted{int(use_cor)}__diff_{srges_diff_mode}"
    _print_metrics(
        f"{disease} | full sRGES | use_cor={use_cor} | mode={srges_diff_mode}",
        metrics,
    )
    _maybe_save(merged, metrics, out_dir=out_dir, disease=disease, tag=tag,
                title=f"{disease} full sRGES  (use_cor={use_cor})")
    return merged, metrics


def run_full_sRGES_netcos(disease, cs_files=None, *,
                          landmark_disease_pref=False,
                          use_cor=False,
                          srges_diff_mode=conf.srges_diff_mode,
                          out_dir=None):
    """Scope C: NetCoS RGES -> my sRGES vs BinChen sRGES (end-to-end test).

    Discrepancies here include both per-signature RGES drift and the sRGES
    aggregation drift (NaN handling etc.).
    """
    if cs_files is None:
        cs_filename = pick_canonical_cs_run(disease, landmark_disease_pref)
        if cs_filename is None:
            raise FileNotFoundError(
                f"No canonical NetCoS CS run found for {disease}. "
                f"Either pass --netcos-cs-file explicitly, or run cs_batch.py "
                f"with the canonical config first."
            )
        cs_files = [cs_out_dir_for(disease, landmark_disease_pref) / cs_filename]
    cs_files = [Path(p) for p in (cs_files if isinstance(cs_files, list) else [cs_files])]

    print(f"\n[run_full_sRGES_netcos] disease={disease}, CS files:")
    for p in cs_files:
        print(f"    {p}")

    inp = (load_netcos_cs_combined(cs_files) if len(cs_files) > 1
           else load_netcos_cs(cs_files[0]))

    ref = load_reference_full_sRGES(disease)

    weights = cor_weights_for_disease(disease) if use_cor else None
    my = compute_sRGES(
        inp, score_col="connectivity_score",
        cell_lines=None,
        cell_line_weights=weights,
        apply_cor_weight=use_cor,
        filter_is_gold=False,
        srges_diff_mode=srges_diff_mode,
    )

    merged, metrics = compare(my, ref, key_col="pert_iname")

    tag = f"netcos_full__corweighted{int(use_cor)}__diff_{srges_diff_mode}"
    _print_metrics(
        f"{disease} | NetCoS-input full sRGES | use_cor={use_cor} | "
        f"mode={srges_diff_mode}",
        metrics,
    )
    _maybe_save(
        merged, metrics, out_dir=out_dir, disease=disease, tag=tag,
        title=f"{disease} NetCoS->sRGES vs BinChen sRGES  (use_cor={use_cor})",
    )
    return merged, metrics


def run_fig3_sRGES(disease, *, cell_line=None, use_cor=False,
                   srges_diff_mode=conf.srges_diff_mode, out_dir=None):
    """Scope B: BinChen RGES restricted to one cell line (Fig.3)."""
    cell_line = cell_line or conf.DEFAULT_CELL_LINE[disease]
    inp = load_binchen_all_lincs_score(disease)
    ref = load_reference_fig3_sRGES(disease)

    weights = cor_weights_for_disease(disease) if use_cor else None
    my = compute_sRGES(
        inp, score_col="RGES",
        cell_lines=[cell_line],
        cell_line_weights=weights,
        apply_cor_weight=use_cor,
        filter_is_gold=False,
        srges_diff_mode=srges_diff_mode,
    )

    merged, metrics = compare(my, ref, key_col="pert_iname")

    tag = f"fig3_{cell_line}__corweighted{int(use_cor)}__diff_{srges_diff_mode}"
    _print_metrics(
        f"{disease} | per-cell-line ({cell_line}) sRGES | "
        f"use_cor={use_cor} | mode={srges_diff_mode}",
        metrics,
    )
    _maybe_save(
        merged, metrics, out_dir=out_dir, disease=disease, tag=tag,
        title=f"{disease} per-cell-line sRGES ({cell_line}, use_cor={use_cor})",
    )
    return merged, metrics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--disease",
                        choices=list(conf.DEFAULT_CELL_LINE.keys()),
                        default="LIHC")
    parser.add_argument("--scope",
                        choices=["full", "per_cell_line", "full_netcos"],
                        default="full",
                        help="full = BinChen RGES -> my sRGES vs BinChen sRGES "
                             "(tests aggregation only). "
                             "per_cell_line = same but restricted to canonical "
                             "single cell line (Fig.3). "
                             "full_netcos = NetCoS RGES -> my sRGES vs BinChen "
                             "sRGES (end-to-end test).")
    parser.add_argument("--cell-line", default=None,
                        help="Cell line for the per_cell_line scope "
                             "(default per-disease canonical from conf.DEFAULT_CELL_LINE).")
    parser.add_argument("--netcos-cs-file", action="append", default=None,
                        help="(scope=full_netcos) Path to a NetCoS CS .tsv file. "
                             "Repeat to combine per-cell-line CS files. If omitted, "
                             "the canonical Bin-Chen-config run is auto-picked from "
                             "logs/cs_runs.tsv.")
    parser.add_argument("--landmark-disease-pref", action="store_true",
                        help="(scope=full_netcos auto-pick) Prefer the *_LM disease run.")
    parser.add_argument("--use-cor", action="store_true",
                        default=conf.srges_use_cor,
                        help="Enable per-cell-line cor*RGES weighting "
                             "(Bin Chen sRGES_all_cmpds.R). "
                             "Default from conf.srges_use_cor.")
    parser.add_argument("--srges-diff-mode",
                        choices=["binchen_r", "corrected"],
                        default=conf.srges_diff_mode)
    parser.add_argument("--out-dir", default=None,
                        help="If set, write joined.tsv / summary.tsv / scatter.png "
                             "there. Defaults to conf.BINCHEN_VALIDATION_OUT_DIR.")
    parser.add_argument("--all", action="store_true",
                        help="Sweep all (disease, scope) combinations.")
    args = parser.parse_args(argv)

    out_dir = (Path(args.out_dir) if args.out_dir
               else conf.BINCHEN_VALIDATION_OUT_DIR)

    if args.all:
        for disease in conf.DEFAULT_CELL_LINE.keys():
            for fn, label in [
                (lambda d: run_full_sRGES(
                    d, use_cor=args.use_cor,
                    srges_diff_mode=args.srges_diff_mode, out_dir=out_dir),
                 "full"),
                (lambda d: run_fig3_sRGES(
                    d, use_cor=args.use_cor,
                    srges_diff_mode=args.srges_diff_mode, out_dir=out_dir),
                 "fig3"),
                (lambda d: run_full_sRGES_netcos(
                    d, landmark_disease_pref=args.landmark_disease_pref,
                    use_cor=args.use_cor,
                    srges_diff_mode=args.srges_diff_mode, out_dir=out_dir),
                 "full_netcos"),
            ]:
                try:
                    fn(disease)
                except FileNotFoundError as e:
                    print(f"  [skip {label} {disease}] missing data: {e}")
        return

    if args.scope == "full":
        run_full_sRGES(args.disease, use_cor=args.use_cor,
                       srges_diff_mode=args.srges_diff_mode, out_dir=out_dir)
    elif args.scope == "full_netcos":
        run_full_sRGES_netcos(args.disease,
                              cs_files=args.netcos_cs_file,
                              landmark_disease_pref=args.landmark_disease_pref,
                              use_cor=args.use_cor,
                              srges_diff_mode=args.srges_diff_mode,
                              out_dir=out_dir)
    else:
        run_fig3_sRGES(args.disease, cell_line=args.cell_line,
                       use_cor=args.use_cor,
                       srges_diff_mode=args.srges_diff_mode,
                       out_dir=out_dir)


if __name__ == "__main__":
    main()
