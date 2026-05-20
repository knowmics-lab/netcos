#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026-05-11

@author: los4

Per-signature RGES replication check against Bin Chen 2017.

For one disease at a time, compare a NetCoS connectivity-score run against
Bin Chen 2017's reference per-LINCS-signature RGES values
(`data/BinChen2017/data/data/<DISEASE>/lincs_score_0.csv` or `_1.csv`).

Join key:
    NetCoS  output column `LINCS_id` (integer LINCS signature id)
        <=>
    BinChen lincs_score_0.csv column `id`

This is the "Figure 3 level" comparison: plain RGES per signature, before any
per-drug sRGES collapse. The corresponding canonical NetCoS hyperparameter set
(see Step 1 notes) is:
    cs_mith=0, CS_ON_PATHWAYS=False, CS_METHOD='bin_chen',
    cs_on_LM=1, landmark_drug=True
which produces a CS file named like '*_DEG_LM_bin_chen_connectivity_score.tsv'.

Outputs are per (disease, cs_run_id, ref_file, cell-line scope):
    - <OOUT_DIR>/binchen2017_validation/<disease>__<cs_run_id>__vs__<ref>__<scope>__summary.tsv
    - <IMG_DIR>/binchen2017_validation/<disease>__<cs_run_id>__vs__<ref>__<scope>__scatter.png
    - <OOUT_DIR>/binchen2017_validation/<disease>__<cs_run_id>__vs__<ref>__<scope>__joined.tsv

Usage:
    # Compare the canonical LIHC DEG-LM-bin_chen run against lincs_score_0 (HEPG2 only):
    python -m src.validations.binchen2017.compare_rges_per_disease \
        --disease LIHC \
        --cs-file '08_05_2026_10_26_DEG_LM_bin_chen_connectivity_score.tsv' \
        --cell-line HEPG2

    # Or auto-pick the canonical run for LIHC from cs_runs.tsv:
    python -m src.validations.binchen2017.compare_rges_per_disease --disease LIHC --auto

    # Compare against the more-selective lincs_score_1.csv reference:
    python -m src.validations.binchen2017.compare_rges_per_disease \
        --disease LIHC --auto --ref lincs_score_1

    # Don't filter by cell_id (compare against ALL LINCS cell lines in the reference):
    python -m src.validations.binchen2017.compare_rges_per_disease \
        --disease LIHC --auto --cell-line ALL
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
from loader import cs_out_dir_for, pick_canonical_cs_run as _pick_canonical_cs_run_shared  # noqa: E402


# ---------------------------------------------------------------------------
# Per-disease cell line mapping for Bin Chen 2017 Figure 3
# ---------------------------------------------------------------------------
DEFAULT_CELL_LINE = {
    "LIHC": "HEPG2",
    "BRCA": "MCF7",
    "COAD": "HT29",
    "ER": "MCF7",   # ER not in SD3/SD4; default to MCF7 (shared lineage with BRCA)
}

# Canonical NetCoS RGES-replication signature inside the cs_run_id filename.
# The make_cs_filename() helper builds names like:
#     <timestamp>_<mith|DEG>[_LM]_[pw_]<CS_METHOD>_connectivity_score.tsv
# So the canonical Bin-Chen RGES replication run contains '_DEG_LM_bin_chen_'
# and crucially NOT 'disease_sorted'.
CANONICAL_TAG = "DEG_LM_bin_chen"


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
# cs_out_dir_for() and pick_canonical_cs_run() now live in src/loader.py so
# this script and test_sRGES_replication.py share one implementation. The
# imports are at the top of the file. We keep a thin local wrapper around
# pick_canonical_cs_run that adds a fallback: if cs_runs.tsv lookup yields
# nothing, glob the local CS output dir for files containing CANONICAL_TAG
# (e.g. '_DEG_LM_bin_chen_'). This fallback predates cs_runs.tsv and
# is kept for backward compat.

def binchen_disease_dir(disease):
    """Return Bin Chen's per-disease data dir on local disk."""
    return conf.BC_DISEASE_DATA / disease


def pick_canonical_cs_run(disease, landmark_disease_pref=False):
    """Pick the most recent NetCoS CS .tsv matching the Bin Chen canonical
    config for `disease`. Delegates to loader.pick_canonical_cs_run, then
    falls back to a CANONICAL_TAG glob if cs_runs.tsv lookup returns None.
    """
    chosen = _pick_canonical_cs_run_shared(
        disease, landmark_disease_pref=landmark_disease_pref
    )
    if chosen is not None:
        return chosen
    # Fallback: glob the local CS output dir for canonical-tag files.
    candidates = sorted(
        cs_out_dir_for(disease, landmark_disease_pref).glob(
            f"*_{CANONICAL_TAG}_connectivity_score.tsv"
        )
    )
    # Filter out 'disease_sorted' which shares the 'bin_chen' substring
    candidates = [c for c in candidates if "disease_sorted" not in c.name]
    if not candidates:
        return None
    return candidates[-1].name


def load_netcos_cs(cs_file_path):
    """Load a NetCoS connectivity_score.tsv produced by cs_batch.py."""
    df = pd.read_csv(cs_file_path, sep="\t", dtype={"LINCS_id": str, "pert_id": str})
    required = {"LINCS_id", "connectivity_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"NetCoS CS file missing columns {missing}: {cs_file_path}")
    return df


def load_binchen_ref(disease, ref_name="lincs_score_0"):
    """
    Load Bin Chen's per-signature RGES reference for `disease`.

    ref_name: 'lincs_score_0' (all signatures) or 'lincs_score_1' (filtered subset).
    """
    path = binchen_disease_dir(disease) / f"{ref_name}.csv"
    df = pd.read_csv(path, dtype={"id": str, "pert_id": str})
    # BinChen CSVs have an unnamed row-index column quoted as "" — pandas keeps it
    # as "Unnamed: 0". Drop it if present.
    if df.columns[0].startswith("Unnamed"):
        df = df.drop(columns=df.columns[0])
    return df


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def topk_overlap(my_series, ref_series, k, ascending=True):
    """Fraction of the top-k items in each ranking that appear in the other."""
    my_top = set(my_series.sort_values(ascending=ascending).head(k).index)
    ref_top = set(ref_series.sort_values(ascending=ascending).head(k).index)
    if not my_top or not ref_top:
        return float("nan")
    return len(my_top & ref_top) / k


def compute_metrics(joined, my_col, ref_col):
    """Compute Pearson, Spearman, Kendall + N for two columns of `joined`."""
    sub = joined[[my_col, ref_col]].dropna()
    out = {
        "n_compared": int(len(sub)),
        "n_my_with_value": int(joined[my_col].notna().sum()),
        "n_ref_with_value": int(joined[ref_col].notna().sum()),
    }
    if len(sub) < 3:
        out.update({"pearson_r": np.nan, "pearson_p": np.nan,
                    "spearman_r": np.nan, "spearman_p": np.nan})
        return out
    pr, pp = stats.pearsonr(sub[my_col], sub[ref_col])
    sr, sp = stats.spearmanr(sub[my_col], sub[ref_col])
    out.update({"pearson_r": pr, "pearson_p": pp,
                "spearman_r": sr, "spearman_p": sp})
    return out


def make_scatter(joined, my_col, ref_col, title, out_png):
    sub = joined[[my_col, ref_col]].dropna()
    pr = sp = np.nan
    if len(sub) >= 3:
        pr, _ = stats.pearsonr(sub[my_col], sub[ref_col])
        sp, _ = stats.spearmanr(sub[my_col], sub[ref_col])

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(sub[ref_col], sub[my_col], s=6, alpha=0.4, edgecolors="none")
    # y=x reference
    lo = float(np.nanmin([sub[my_col].min(), sub[ref_col].min()]))
    hi = float(np.nanmax([sub[my_col].max(), sub[ref_col].max()]))
    ax.plot([lo, hi], [lo, hi], "--", color="gray", linewidth=0.8, label="y = x")
    # OLS fit line
    if len(sub) >= 2:
        slope, intercept, _, _, _ = stats.linregress(sub[ref_col], sub[my_col])
        xs = np.array([lo, hi])
        ax.plot(xs, slope * xs + intercept, "-", color="C3", linewidth=1.0,
                label=f"OLS: y = {slope:.3f} x + {intercept:.3f}")
    ax.set_xlabel(f"Bin Chen 2017 — {ref_col}")
    ax.set_ylabel(f"NetCoS — {my_col}")
    ax.set_title(
        f"{title}\nN={len(sub)}  Pearson r={pr:.3f}  Spearman ρ={sp:.3f}"
    )
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Core comparison
# ---------------------------------------------------------------------------
def compare(disease, cs_file_name, ref_name="lincs_score_0", cell_line=None,
            landmark_disease=False, topk=(10, 50, 100, 500),
            img_dir=None, oout_dir=None, verbose=True):
    """
    Run one comparison and write summary, joined table, and scatter plot.

    Returns the summary dict.
    """
    # --- paths ---
    cs_out_d = cs_out_dir_for(disease, landmark_disease=landmark_disease)
    cs_path = cs_out_d / cs_file_name
    if not cs_path.exists():
        raise FileNotFoundError(f"NetCoS CS output not found: {cs_path}")

    cs_run_id = Path(cs_file_name).stem  # drop .tsv

    if img_dir is None:
        img_dir = conf.IMG_DIR / "binchen2017_validation"
    if oout_dir is None:
        oout_dir = conf.OOUT_DIR / "binchen2017_validation"
    img_dir.mkdir(parents=True, exist_ok=True)
    oout_dir.mkdir(parents=True, exist_ok=True)

    scope_tag = (cell_line or "ALL").upper()
    base_name = f"{disease}__{cs_run_id}__vs__{ref_name}__{scope_tag}"
    summary_path = oout_dir / f"{base_name}__summary.tsv"
    joined_path  = oout_dir / f"{base_name}__joined.tsv"
    scatter_path = img_dir  / f"{base_name}__scatter.png"

    # --- load ---
    if verbose:
        print(f"[load] NetCoS CS: {cs_path}")
    mine = load_netcos_cs(cs_path)
    if verbose:
        print(f"       {len(mine):,} rows, {mine['LINCS_id'].nunique():,} unique LINCS_id")

    if verbose:
        print(f"[load] Bin Chen ref: {binchen_disease_dir(disease) / (ref_name + '.csv')}")
    ref = load_binchen_ref(disease, ref_name=ref_name)
    if verbose:
        print(f"       {len(ref):,} rows, {ref['id'].nunique():,} unique id; "
              f"cell_ids: {sorted(ref['cell_id'].dropna().unique())[:8]}"
              f"{' ...' if ref['cell_id'].nunique() > 8 else ''}")

    # --- optional cell-line filter on Bin Chen side ---
    if cell_line and cell_line.upper() != "ALL":
        n_before = len(ref)
        ref = ref[ref["cell_id"].astype(str).str.upper() == cell_line.upper()].copy()
        if verbose:
            print(f"[filter] cell_id == {cell_line.upper()}: "
                  f"{n_before:,} -> {len(ref):,} rows")

    # --- join ---
    # Coerce both keys to str to avoid float<->int mismatches.
    mine_j = mine.rename(columns={"LINCS_id": "_key"}).copy()
    ref_j  = ref.rename(columns={"id": "_key"}).copy()
    mine_j["_key"] = mine_j["_key"].astype(str)
    ref_j["_key"]  = ref_j["_key"].astype(str)

    # Prefix every Bin Chen column with 'ref_' (except the join key) so
    # downstream code can reference them as ref_RGES, ref_pearson, ref_cell_id,
    # ref_pert_iname, etc. regardless of whether the column name also exists
    # in the NetCoS output.
    ref_j = ref_j.rename(
        columns={c: f"ref_{c}" for c in ref_j.columns if c != "_key"}
    )

    joined = mine_j.merge(ref_j, on="_key", how="inner")
    if verbose:
        n_mine = len(mine_j)
        n_ref  = len(ref_j)
        print(f"[join] inner-join on LINCS_id: my N={n_mine:,}, ref N={n_ref:,}, "
              f"joined N={len(joined):,} "
              f"(unmatched_my={n_mine - joined['_key'].nunique():,}, "
              f"unmatched_ref={n_ref - joined['_key'].nunique():,})")

    if len(joined) == 0:
        print("[error] No rows matched after join. Aborting comparison.")
        return None

    # --- metrics ---
    # NetCoS connectivity_score should align with BinChen RGES (Pearson/Spearman).
    metrics = {}
    metrics["cs_vs_rges"] = compute_metrics(joined, "connectivity_score", "ref_RGES")
    # Sanity-check Pearson/Spearman correlations of underlying signatures
    # (the CS file also stores pearson/spearman per signature; BinChen's ref has the same).
    if "pearson" in joined.columns and "ref_pearson" in joined.columns:
        metrics["my_pearson_vs_ref_pearson"] = compute_metrics(
            joined, "pearson", "ref_pearson"
        )
    if "spearman" in joined.columns and "ref_spearman" in joined.columns:
        metrics["my_spearman_vs_ref_spearman"] = compute_metrics(
            joined, "spearman", "ref_spearman"
        )

    # Top/bottom-K overlap (rank lowest = most reversed signatures)
    overlap_rows = []
    joined_for_topk = joined.set_index("_key")
    for k in topk:
        if len(joined_for_topk) < k:
            continue
        ov_low = topk_overlap(
            joined_for_topk["connectivity_score"], joined_for_topk["ref_RGES"],
            k=k, ascending=True
        )
        ov_hi = topk_overlap(
            joined_for_topk["connectivity_score"], joined_for_topk["ref_RGES"],
            k=k, ascending=False
        )
        overlap_rows.append({"k": k, "topk_low_overlap": ov_low, "topk_high_overlap": ov_hi})

    # --- write summary ---
    summary_rows = []
    for label, m in metrics.items():
        summary_rows.append({
            "comparison": label,
            **m,
        })
    summary_df = pd.DataFrame(summary_rows)
    overlap_df = pd.DataFrame(overlap_rows)
    overlap_df.insert(0, "comparison", "topk_overlap(cs_vs_rges)")
    # Concatenate into one TSV with a comparison column distinguishing them.
    out_tbl = pd.concat([summary_df, overlap_df], axis=0, ignore_index=True, sort=False)

    # Stamp the run with provenance.
    header_rows = pd.DataFrame([
        {"comparison": "_meta:disease",      "value": disease},
        {"comparison": "_meta:cs_run_id",    "value": cs_run_id},
        {"comparison": "_meta:ref_file",     "value": ref_name},
        {"comparison": "_meta:cell_scope",   "value": scope_tag},
        {"comparison": "_meta:n_my_rows",    "value": len(mine_j)},
        {"comparison": "_meta:n_ref_rows",   "value": len(ref_j)},
        {"comparison": "_meta:n_joined",     "value": len(joined)},
    ])
    out_tbl = pd.concat([header_rows, out_tbl], axis=0, ignore_index=True, sort=False)
    out_tbl.to_csv(summary_path, sep="\t", index=False)

    # --- write joined table (compact, useful for follow-up analyses) ---
    keep_cols = ["_key", "connectivity_score", "ref_RGES",
                 "pearson", "ref_pearson", "spearman", "ref_spearman",
                 "pert_id", "ref_pert_id", "ref_pert_iname",
                 "ref_cell_id", "ref_pert_dose", "ref_pert_time", "ref_is_gold"]
    keep_cols = [c for c in keep_cols if c in joined.columns]
    joined[keep_cols].rename(columns={"_key": "LINCS_id"}).to_csv(
        joined_path, sep="\t", index=False
    )

    # --- plot ---
    make_scatter(
        joined,
        my_col="connectivity_score",
        ref_col="ref_RGES",
        title=f"{disease} ({scope_tag}) — NetCoS CS vs Bin Chen RGES\n"
              f"run: {cs_run_id}   ref: {ref_name}",
        out_png=scatter_path,
    )

    if verbose:
        print()
        print(f"[result] CS vs RGES: N={metrics['cs_vs_rges']['n_compared']:,}  "
              f"Pearson r={metrics['cs_vs_rges']['pearson_r']:.4f}  "
              f"Spearman ρ={metrics['cs_vs_rges']['spearman_r']:.4f}")
        if "my_pearson_vs_ref_pearson" in metrics:
            print(f"         my-pearson vs ref-pearson Pearson r="
                  f"{metrics['my_pearson_vs_ref_pearson']['pearson_r']:.4f} "
                  f"(should be ~1 for a faithful replication)")
        for row in overlap_rows:
            print(f"         top-{row['k']:<4} low overlap (most reversed): "
                  f"{row['topk_low_overlap']:.3f}")
        print(f"[write] summary   -> {summary_path}")
        print(f"[write] joined    -> {joined_path}")
        print(f"[write] scatter   -> {scatter_path}")

    return {
        "disease": disease,
        "cs_run_id": cs_run_id,
        "ref": ref_name,
        "cell_scope": scope_tag,
        "metrics": metrics,
        "topk": overlap_rows,
        "summary_path": str(summary_path),
        "joined_path": str(joined_path),
        "scatter_path": str(scatter_path),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compare a NetCoS per-signature CS run against Bin Chen 2017's "
            "per-signature RGES reference (lincs_score_0.csv or lincs_score_1.csv) "
            "for one disease."
        )
    )
    p.add_argument("--disease", default=conf.DISEASE,
                   help="Disease symbol (LIHC/BRCA/COAD/ER). "
                        f"Default from conf.py: {conf.DISEASE}")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--cs-file", default=None,
                   help="NetCoS CS output filename (just the basename, located in "
                        "connectivity_score/output/<DISEASE>{,_LM}/). "
                        "If omitted, use --auto.")
    g.add_argument("--auto", action="store_true",
                   help="Auto-pick the most recent canonical Bin Chen RGES run "
                        "(mith=0, cs_on_LM=1, CS_ON_PATHWAYS=0, CS_METHOD=bin_chen) "
                        "from logs/cs_runs.tsv.")
    p.add_argument("--ref", default="lincs_score_1",
                   choices=["lincs_score_0", "lincs_score_1"],
                   help="Bin Chen reference file (per-disease folder). "
                        "lincs_score_0 = all signatures; lincs_score_1 = filtered subset.")
    p.add_argument("--cell-line", default=None,
                   help="Filter Bin Chen reference by cell_id. Default: per-disease "
                        "canonical (LIHC=HEPG2, BRCA=MCF7, COAD=HT29, ER=MCF7). "
                        "Pass 'ALL' to skip filtering.")
    p.add_argument("--landmark-disease", action="store_true",
                   help="Look for the CS run under the *_LM disease output directory "
                        "(set this if landmark_disease=True in the run).")
    return p.parse_args()


def main():
    args = _parse_args()
    disease = args.disease

    # Resolve which CS file to compare.
    cs_file = args.cs_file
    if cs_file is None:
        cs_file = pick_canonical_cs_run(disease, landmark_disease_pref=args.landmark_disease)
        if cs_file is None:
            print(f"[error] No canonical Bin Chen RGES CS run found for {disease}. "
                  f"Provide --cs-file explicitly.")
            sys.exit(2)
        print(f"[auto-pick] {cs_file}")

    # Resolve cell-line scope.
    cell_line = args.cell_line
    if cell_line is None:
        cell_line = DEFAULT_CELL_LINE.get(disease)
        if cell_line is None:
            print(f"[warn] No default cell line known for {disease}; falling back to ALL.")
            cell_line = "ALL"

    compare(
        disease=disease,
        cs_file_name=cs_file,
        ref_name=args.ref,
        cell_line=cell_line,
        landmark_disease=args.landmark_disease,
    )


if __name__ == "__main__":
    main()
