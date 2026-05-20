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

The cor-weighted variant (sRGES_all_cmpds.R) can be enabled via --use-cor.

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
from sRGES import compute_sRGES, load_cell_line_weights_binchen  # noqa: E402
from loader import cs_out_dir_for, pick_canonical_cs_run  # noqa: E402


DEFAULT_CELL_LINE = {
    "LIHC": "HEPG2",
    "BRCA": "MCF7",
    "COAD": "HT29",
    "ER":   "MCF7",
}


def binchen_disease_dir(disease):
    return conf.BC_DISEASE_DATA / disease


# cs_out_dir_for and pick_canonical_cs_run are imported from src/loader.py
# (single source of truth shared with compare_rges_per_disease.py).


def add_lincs_metadata_to_cs(
    cs_df,
    lincs_metadata_path,
    cs_id_col="LINCS_id",
    metadata_id_col="id",
    extra_cols=("pert_iname", "cell_id", "pert_dose", "pert_time", "is_gold"),
):
    """Enrich a NetCoS CS dataframe with the LINCS metadata columns needed by sRGES.

    Joins on LINCS_id <-> id in lincs_sig_info{_new}.csv.
    """
    usecols = [metadata_id_col, *extra_cols]
    meta = pd.read_csv(lincs_metadata_path, usecols=usecols, dtype="str") \
             .drop_duplicates(subset=[metadata_id_col])
    out = cs_df.merge(meta, left_on=cs_id_col, right_on=metadata_id_col, how="left")
    if metadata_id_col in out.columns and metadata_id_col != cs_id_col:
        out = out.drop(columns=[metadata_id_col])
    # Coerce dose/time to numeric so compute_sRGES filters work.
    if "pert_dose" in out.columns:
        out["pert_dose"] = pd.to_numeric(out["pert_dose"], errors="coerce")
    if "pert_time" in out.columns:
        out["pert_time"] = pd.to_numeric(out["pert_time"], errors="coerce")
    n_missing = out["pert_iname"].isna().sum() if "pert_iname" in out.columns else 0
    if n_missing:
        print(f"[warn] add_lincs_metadata_to_cs: {n_missing} LINCS_id rows had "
              f"no metadata match (dropped from sRGES input).")
        out = out.dropna(subset=["pert_iname", "cell_id", "pert_dose", "pert_time"])
    return out


def load_netcos_cs(cs_file_path, lincs_metadata_path):
    """Load a NetCoS CS .tsv and enrich with LINCS metadata (per-signature)."""
    df = pd.read_csv(cs_file_path, sep="\t", dtype={"LINCS_id": str})
    if "LINCS_id" not in df.columns or "connectivity_score" not in df.columns:
        raise ValueError(
            f"NetCoS CS file missing required columns LINCS_id/connectivity_score: "
            f"{cs_file_path}"
        )
    return add_lincs_metadata_to_cs(df, lincs_metadata_path)


def load_netcos_cs_combined(cs_file_paths, lincs_metadata_path):
    """Concatenate multiple per-cell-line NetCoS CS .tsv files and enrich.

    Useful when the user's pipeline produces one CS file per cell-line drug
    set (HEPG2 + MCF7 + HT29 + ...). The combined dataframe is what BinChen's
    full sRGES expects (one row per (drug-signature, cell)).
    """
    parts = []
    for p in cs_file_paths:
        d = pd.read_csv(p, sep="\t", dtype={"LINCS_id": str})
        d["__cs_source_file__"] = Path(p).name
        parts.append(d)
    if not parts:
        raise ValueError("No CS files provided")
    df = pd.concat(parts, ignore_index=True)
    # Same LINCS_id can appear in multiple per-cell-line files only if BinChen
    # ever cross-scored a signature against multiple disease cell lines;
    # drop_duplicates(subset='LINCS_id') is intentionally NOT applied because
    # for sRGES we want every (sig, cell) combination preserved.
    return add_lincs_metadata_to_cs(df, lincs_metadata_path)


# _parse_cs_run_id_ts and pick_canonical_cs_run live in src/loader.py and
# are imported at the top of this module. Shared with
# compare_rges_per_disease.py so both replications use the same picker.


def load_binchen_all_lincs_score(disease):
    """Load Bin Chen's per-instance RGES file (their per-signature output).

    Note: the column is named `cmap_score` in this file, but semantically it
    is the per-instance RGES (matches `RGES` in lincs_score_1.csv exactly for
    rows that overlap).
    """
    path = binchen_disease_dir(disease) / "all_lincs_score.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    df = df.rename(columns={"cmap_score": "RGES"})
    # Coerce dose/time to numeric
    df["pert_dose"] = pd.to_numeric(df["pert_dose"], errors="coerce")
    df["pert_time"] = pd.to_numeric(df["pert_time"], errors="coerce")
    return df


def load_reference_full_sRGES(disease):
    """Load Bin Chen's per-drug sRGES (SD6-equivalent local file)."""
    path = binchen_disease_dir(disease) / "lincs_cancer_sRGES.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    return df  # cols: index, pert_iname, mean, n, median, sd, sRGES


def load_reference_fig3_sRGES(disease):
    """Load Bin Chen's per-drug sRGES restricted to drugs with IC50 (SD5 / Fig.3)."""
    path = binchen_disease_dir(disease) / "rges_ic50_normalized.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    return df  # cols: index, pert_iname, sRGES, standard_value, activity


def cor_weights_for_disease(disease):
    """Build the lincs_cell_id -> cor map a la sRGES_all_cmpds.R."""
    cl_path = binchen_disease_dir(disease) / f"cell_line_{disease}_tacle.csv"
    ccle_path = conf.BC_DISEASE_DATA / "raw" / "cell_line_lincs_ccle.csv"
    if not cl_path.exists() or not ccle_path.exists():
        return None
    return load_cell_line_weights_binchen(str(cl_path), str(ccle_path))


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
    mine = merged[my_score_col + "_mine"] if my_score_col + "_mine" in merged.columns else merged[my_score_col]
    refv = merged[ref_score_col + "_ref"] if ref_score_col + "_ref" in merged.columns else merged[ref_score_col]

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


# ---------------------------------------------------------------------------
# Main per-disease drivers
# ---------------------------------------------------------------------------

def run_full_sRGES(disease, *, use_cor=False, srges_diff_mode="binchen_r",
                   out_dir=None):
    inp = load_binchen_all_lincs_score(disease)
    ref = load_reference_full_sRGES(disease)

    weights = cor_weights_for_disease(disease) if use_cor else None
    my = compute_sRGES(
        inp, score_col="RGES", drug_col="pert_iname",
        cell_col="cell_id", dose_col="pert_dose", time_col="pert_time",
        is_gold_col="is_gold", filter_is_gold=False,
        cell_lines=None,  # full -> all cell lines
        cell_line_weights=weights,
        apply_cor_weight=use_cor,
        srges_diff_mode=srges_diff_mode,
    )

    merged, metrics = compare(my, ref, key_col="pert_iname")

    tag = f"full__corweighted{int(use_cor)}__diff_{srges_diff_mode}"
    print(f"\n=== {disease} | full sRGES | use_cor={use_cor} | "
          f"mode={srges_diff_mode} ===")
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

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out_dir / f"{disease}__{tag}__joined.tsv",
                      sep="\t", index=False)
        pd.DataFrame([metrics]).to_csv(
            out_dir / f"{disease}__{tag}__summary.tsv", sep="\t", index=False
        )
        plot_scatter(
            merged, my_col="sRGES_mine", ref_col="sRGES_ref",
            title=f"{disease} full sRGES  (use_cor={use_cor})",
            out_path=out_dir / f"{disease}__{tag}__scatter.png",
        )
    return merged, metrics


def run_full_sRGES_netcos(disease, cs_files=None, *,
                           landmark_disease_pref=False,
                           use_cor=False, srges_diff_mode="binchen_r",
                           out_dir=None):
    """End-to-end NetCoS test: run my compute_sRGES on a NetCoS CS file (or list
    of per-cell-line CS files), and compare to BinChen's lincs_cancer_sRGES.csv.

    This is the stronger replication: NetCoS RGES -> my sRGES vs BinChen RGES
    -> BinChen sRGES. Discrepancies here include both the per-signature RGES
    drift (which is bit-exact modulo the sign-flip fix) and the sRGES
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

    inp = (load_netcos_cs_combined(cs_files, conf.lincs_metadata_path)
           if len(cs_files) > 1
           else load_netcos_cs(cs_files[0], conf.lincs_metadata_path))

    ref = load_reference_full_sRGES(disease)

    weights = cor_weights_for_disease(disease) if use_cor else None
    my = compute_sRGES(
        inp, score_col="connectivity_score", drug_col="pert_iname",
        cell_col="cell_id", dose_col="pert_dose", time_col="pert_time",
        is_gold_col="is_gold", filter_is_gold=False,
        cell_lines=None,
        cell_line_weights=weights,
        apply_cor_weight=use_cor,
        srges_diff_mode=srges_diff_mode,
    )

    merged, metrics = compare(my, ref, key_col="pert_iname")

    tag = (f"netcos_full__corweighted{int(use_cor)}__"
           f"diff_{srges_diff_mode}")
    print(f"\n=== {disease} | NetCoS-input full sRGES | use_cor={use_cor} | "
          f"mode={srges_diff_mode} ===")
    print(f"  n joined drugs                 : {metrics['n_joined']}")
    print(f"  Pearson  r vs BinChen sRGES    : {metrics['pearson_r']:.6f}  "
          f"(p={metrics['pearson_p']:.2e})")
    print(f"  Spearman rho vs BinChen sRGES  : {metrics['spearman_rho']:.6f}  "
          f"(p={metrics['spearman_p']:.2e})")
    print(f"  RMSE / MAE                     : {metrics['rmse']:.6f} / "
          f"{metrics['mae']:.6f}")
    for k in (10, 50, 100):
        if f"top{k}_overlap" in metrics:
            print(f"  top-{k} overlap (low)           : {metrics[f'top{k}_overlap']:.3f}")

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out_dir / f"{disease}__{tag}__joined.tsv",
                      sep="\t", index=False)
        pd.DataFrame([metrics]).to_csv(
            out_dir / f"{disease}__{tag}__summary.tsv", sep="\t", index=False
        )
        plot_scatter(
            merged, my_col="sRGES_mine", ref_col="sRGES_ref",
            title=f"{disease} NetCoS->sRGES vs BinChen sRGES  (use_cor={use_cor})",
            out_path=out_dir / f"{disease}__{tag}__scatter.png",
        )
    return merged, metrics


def run_fig3_sRGES(disease, *, cell_line=None, use_cor=False,
                   srges_diff_mode="binchen_r", out_dir=None):
    cell_line = cell_line or DEFAULT_CELL_LINE[disease]
    inp = load_binchen_all_lincs_score(disease)
    ref = load_reference_fig3_sRGES(disease)

    weights = cor_weights_for_disease(disease) if use_cor else None
    my = compute_sRGES(
        inp, score_col="RGES", drug_col="pert_iname",
        cell_col="cell_id", dose_col="pert_dose", time_col="pert_time",
        is_gold_col="is_gold", filter_is_gold=False,
        cell_lines=[cell_line],
        cell_line_weights=weights,
        apply_cor_weight=use_cor,
        srges_diff_mode=srges_diff_mode,
    )

    merged, metrics = compare(my, ref, key_col="pert_iname")

    tag = (f"fig3_{cell_line}__corweighted{int(use_cor)}__"
           f"diff_{srges_diff_mode}")
    print(f"\n=== {disease} | per-cell-line ({cell_line}) sRGES | "
          f"use_cor={use_cor} | mode={srges_diff_mode} ===")
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

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out_dir / f"{disease}__{tag}__joined.tsv",
                      sep="\t", index=False)
        pd.DataFrame([metrics]).to_csv(
            out_dir / f"{disease}__{tag}__summary.tsv", sep="\t", index=False
        )
        plot_scatter(
            merged, my_col="sRGES_mine", ref_col="sRGES_ref",
            title=f"{disease} per-cell-line sRGES ({cell_line}, use_cor={use_cor})",
            out_path=out_dir / f"{disease}__{tag}__scatter.png",
        )
    return merged, metrics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--disease",
                        choices=list(DEFAULT_CELL_LINE.keys()),
                        default="LIHC")
    parser.add_argument("--scope", choices=["full", "per_cell_line", "full_netcos"],
                        default="full",
                        help="full = BinChen RGES -> my sRGES vs BinChen sRGES "
                             "(tests aggregation only). "
                             "per_cell_line = same but restricted to canonical "
                             "single cell line (Fig.3). "
                             "full_netcos = NetCoS RGES -> my sRGES vs BinChen "
                             "sRGES (end-to-end test).")
    parser.add_argument("--cell-line", default=None,
                        help="Cell line for the per_cell_line scope (default per-disease canonical).")
    parser.add_argument("--netcos-cs-file", action="append", default=None,
                        help="(scope=full_netcos) Path to a NetCoS CS .tsv file. "
                             "Repeat to combine per-cell-line CS files. If omitted, "
                             "the canonical Bin-Chen-config run is auto-picked from "
                             "logs/cs_runs.tsv.")
    parser.add_argument("--landmark-disease-pref", action="store_true",
                        help="(scope=full_netcos auto-pick) Prefer the *_LM disease run.")
    parser.add_argument("--use-cor", action="store_true",
                        help="Enable per-cell-line cor*RGES weighting (Bin Chen sRGES_all_cmpds.R).")
    parser.add_argument("--srges-diff-mode",
                        choices=["binchen_r", "corrected"], default="binchen_r")
    parser.add_argument("--out-dir", default=None,
                        help="If set, write joined.tsv / summary.tsv / scatter.png there.")
    parser.add_argument("--all", action="store_true",
                        help="Sweep all (disease, scope) combinations.")
    args = parser.parse_args(argv)

    out_dir = (Path(args.out_dir) if args.out_dir else
               (conf.OTHER_OUT_DIR / "binchen2017_validation"
                if hasattr(conf, "OTHER_OUT_DIR") else None))

    if args.all:
        for disease in DEFAULT_CELL_LINE.keys():
            try:
                run_full_sRGES(disease, use_cor=args.use_cor,
                                srges_diff_mode=args.srges_diff_mode,
                                out_dir=out_dir)
            except FileNotFoundError as e:
                print(f"  [skip full {disease}] missing data: {e}")
            try:
                run_fig3_sRGES(disease, use_cor=args.use_cor,
                                srges_diff_mode=args.srges_diff_mode,
                                out_dir=out_dir)
            except FileNotFoundError as e:
                print(f"  [skip fig3 {disease}] missing data: {e}")
            try:
                run_full_sRGES_netcos(disease,
                                      landmark_disease_pref=args.landmark_disease_pref,
                                      use_cor=args.use_cor,
                                      srges_diff_mode=args.srges_diff_mode,
                                      out_dir=out_dir)
            except FileNotFoundError as e:
                print(f"  [skip full_netcos {disease}] missing data: {e}")
        return

    if args.scope == "full":
        run_full_sRGES(args.disease, use_cor=args.use_cor,
                       srges_diff_mode=args.srges_diff_mode,
                       out_dir=out_dir)
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
