#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sRGES (summarized Reverse Gene Expression Score) — faithful Python port of the
algorithm published by Bin Chen et al. (Nat. Commun. 2017), as implemented in
the Bin-Chen-Lab/RGES public repository:
    https://github.com/Bin-Chen-Lab/RGES/blob/master/sRGES_all_cmpds.R
    https://github.com/Bin-Chen-Lab/RGES/blob/master/runRGESExample.R
    https://github.com/Bin-Chen-Lab/RGES/blob/master/core_functions.R  (getsRGES)

Author: L-F-S, 2026-05-13

Algorithm overview
------------------
For each drug, the per-instance RGES values measured across LINCS L1000
signatures (different cell lines, doses, treatment times) are aggregated into a
single per-drug score, sRGES. The aggregation:

  1. Restricts to per-instance signatures with pert_dose > 0 and
     pert_time in {6, 24} hours (and optionally is_gold == 1).
  2. Builds a dose/time correction model: for each (cell_line, drug) it forms
     all (instance_i, instance_j) pairs where instance_i is a "reference"
     (pert_time == 24h AND pert_dose == 10 µM). The mean RGES difference
     ref - other is computed within each (dose_bin, pert_time) bucket, where
     dose_bin is `low` (<10) or `high` (>=10).
  3. Per-instance adjustment: every individual RGES is shifted by the bucket
     mean diff (which is 0 for "high 24h", i.e. the reference bucket itself).
  4. Optional cell-line weighting: multiply by cor / max(cor) where cor is the
     per-cell-line correlation with the target tumour type (or 1 if unspecified).
  5. Per-drug aggregation: sRGES = mean of adjusted RGES across all of that
     drug's surviving instances.

Notes on faithful reproduction
------------------------------
We reproduce the *exact* indexing pattern used by Bin Chen's R `getsRGES` so
that the Python output matches the published reference outputs (e.g.
data/BinChen2017/data/data/<disease>/lincs_cancer_sRGES.csv).

In the R source, `diff` is a named numeric vector produced by

    tapply(cmap_diff, paste(dose_bin, pert_time), mean)

which R sorts alphabetically as
    diff[1] = "high 24"   (reference; mean ≈ 0)
    diff[2] = "high 6"
    diff[3] = "low 24"
    diff[4] = "low 6"

`getsRGES` then adds:
    short & low  -> diff[4]   ("low 6")   ✓ semantically correct
    long  & low  -> diff[2]   ("high 6")  ← integer-indexed in R source; we mirror
    short & high -> diff[1]   ("high 24") ← integer-indexed in R source; we mirror

We expose `srges_diff_mode` to switch between:
  - 'binchen_r' (default) — bit-exact reproduction of the R behaviour above.
  - 'corrected'           — semantically corrected mapping
                            (long&low -> "low 24", short&high -> "high 6").

NetCoS sign convention
----------------------
NetCoS's `bin_chen` connectivity score now uses the same sign convention as
Bin Chen's RGES (after the sign fix in src/connectivity_score.py::rank_genes,
where the drug signature is sorted descending). Callers can therefore pass
NetCoS CS values directly into compute_sRGES without any sign manipulation.

"""
from __future__ import annotations

import logging
import numpy as np
import pandas as pd
from typing import Iterable, Mapping, Optional, Union
from conf import LINCS_METADATA_PATH

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Step 1: reference dose/time correction model (`diff`)
# ---------------------------------------------------------------------------

def build_reference_diff(
    df: pd.DataFrame,
    *,
    score_col: str = "RGES",
    drug_col: str = "pert_iname",
    cell_col: str = "cell_id",
    dose_col: str = "pert_dose",
    time_col: str = "pert_time",
    ref_dose: float = 10.0,
    ref_time: float = 24.0,
    low_dose_threshold: float = 10.0,
    valid_times: Iterable = (6, 24),
) -> pd.Series:
    """
    Build the (dose_bin, pert_time) -> mean RGES diff lookup table, as in
    sRGES_all_cmpds.R::

        sub <- subset(prediction, pert_dose > 0 & pert_time %in% c(6, 24))
        pairs <- merge(sub, sub, by=c("pert_iname", "cell_id"))
        pairs <- subset(pairs, id.x != id.y &
                                pert_time.x == 24 & pert_dose.x == 10)
        pairs$cmap_diff <- pairs$cmap_score.x - pairs$cmap_score.y
        pairs$dose_bin <- ifelse(pairs$pert_dose.y < 10, "low", "high")
        diff <- tapply(pairs$cmap_diff,
                       paste(pairs$dose_bin, pairs$pert_time.y),
                       mean)

    Returns a pandas Series indexed by the four bucket strings
    "high 24", "high 6", "low 24", "low 6" (alphabetically sorted, matching
    R's tapply-on-character-vector ordering).
    """
    sub = df.copy()
    sub[dose_col] = pd.to_numeric(sub[dose_col], errors="coerce")
    sub[time_col] = pd.to_numeric(sub[time_col], errors="coerce")
    sub = sub[(sub[dose_col] > 0) & (sub[time_col].isin(valid_times))].copy()

    # We need a row-level id to mimic R's `id.x != id.y` filter. The CSV files
    # already have an `id` column (signature id); if it is present we use it,
    # otherwise we fall back to the dataframe index.
    if "id" in sub.columns:
        sub["__row_id__"] = sub["id"]
    else:
        sub["__row_id__"] = sub.index

    # Keep only the columns we need for the self-merge, to keep it light.
    keep = [drug_col, cell_col, "__row_id__", score_col, dose_col, time_col]
    sub_small = sub[keep]

    pairs = sub_small.merge(sub_small, on=[drug_col, cell_col],
                            suffixes=("_x", "_y"))

    # Reference: x side is pert_time == ref_time and pert_dose == ref_dose.
    # Drop trivial self-pairs (id.x == id.y).
    ref_mask = (
        (pairs[f"{dose_col}_x"] == ref_dose)
        & (pairs[f"{time_col}_x"] == ref_time)
        & (pairs["__row_id___x"] != pairs["__row_id___y"])
    )
    pairs = pairs.loc[ref_mask]

    if pairs.empty:
        logger.warning(
            "build_reference_diff: no (drug, cell) pairs found with a "
            "reference instance at dose=%s time=%s. Returning zero diff.",
            ref_dose, ref_time,
        )
        return pd.Series(
            {"high 24": 0.0, "high 6": 0.0, "low 24": 0.0, "low 6": 0.0},
            name="diff",
        )

    cmap_diff = pairs[f"{score_col}_x"].to_numpy() - pairs[f"{score_col}_y"].to_numpy()
    dose_bin = np.where(pairs[f"{dose_col}_y"].to_numpy() < low_dose_threshold,
                        "low", "high")
    # paste(dose_bin, pert_time.y) — R coerces numeric to character without
    # decimals, e.g. 24 -> "24", 6 -> "6"
    time_str = pairs[f"{time_col}_y"].astype(int).astype(str).to_numpy()
    bucket = pd.Series([f"{d} {t}" for d, t in zip(dose_bin, time_str)])

    diff = pd.Series(cmap_diff).groupby(bucket).mean().sort_index()
    diff.name = "diff"
    return diff


# ---------------------------------------------------------------------------
# Step 2: per-instance sRGES adjustment (port of getsRGES)
# ---------------------------------------------------------------------------

def _get_sRGES_one(
    rges: float,
    cor_value: float,
    pert_dose: float,
    pert_time: float,
    diff_lookup: Mapping[str, float],
    max_cor: float,
    *,
    low_dose_threshold: float = 10.0,
    short_time_threshold: float = 24.0,
    apply_cor_weight: bool = True,
    srges_diff_mode: str = "binchen_r",
) -> float:
    """Port of Bin Chen's R::getsRGES (core_functions.R)."""
    if pd.isna(rges) or pd.isna(pert_dose) or pd.isna(pert_time):
        return float("nan")

    out = rges
    bin_time = "short" if pert_time < short_time_threshold else "long"
    bin_dose = "low" if pert_dose < low_dose_threshold else "high"

    if srges_diff_mode == "binchen_r":
        # Match the integer indexing of the R source verbatim:
        #   diff[1]="high 24", diff[2]="high 6", diff[3]="low 24", diff[4]="low 6"
        if bin_time == "short" and bin_dose == "low":
            out += diff_lookup.get("low 6", 0.0)         # diff[4]
        if bin_dose == "low" and bin_time == "long":
            out += diff_lookup.get("high 6", 0.0)        # diff[2]  (per R source)
        if bin_dose == "high" and bin_time == "short":
            out += diff_lookup.get("high 24", 0.0)       # diff[1]  (per R source)
        # high & long: reference, no shift.

    elif srges_diff_mode == "corrected":
        if bin_time == "short" and bin_dose == "low":
            out += diff_lookup.get("low 6", 0.0)
        if bin_dose == "low" and bin_time == "long":
            out += diff_lookup.get("low 24", 0.0)
        if bin_dose == "high" and bin_time == "short":
            out += diff_lookup.get("high 6", 0.0)

    else:
        raise ValueError(
            f"srges_diff_mode must be 'binchen_r' or 'corrected', got "
            f"{srges_diff_mode!r}"
        )

    if apply_cor_weight and max_cor != 0 and not pd.isna(max_cor):
        out = out * (cor_value / max_cor)
    return out


# ---------------------------------------------------------------------------
# Step 3: full pipeline — compute sRGES per drug
# ---------------------------------------------------------------------------

def compute_sRGES(
    rges_df: pd.DataFrame,
    *,
    score_col: str = "RGES",
    drug_col: str = "pert_iname",
    cell_col: str = "cell_id",
    dose_col: str = "pert_dose",
    time_col: str = "pert_time",
    is_gold_col: Optional[str] = "is_gold",
    filter_is_gold: bool = False,
    ref_dose: float = 10.0,
    ref_time: float = 24.0,
    valid_times: Iterable = (6, 24),
    low_dose_threshold: float = 10.0,
    short_time_threshold: float = 24.0,
    cell_lines: Union[None, str, Iterable[str]] = None,
    cell_line_weights: Optional[Mapping[str, float]] = None,
    apply_cor_weight: bool = True,
    srges_diff_mode: str = "binchen_r",
    return_full: bool = True,
) -> pd.DataFrame:
    """
    Compute sRGES per drug.

    Parameters
    ----------
    rges_df : DataFrame
        Per-instance RGES with at least these columns:
            - `score_col` : RGES (or NetCoS connectivity score; see `sign_flip`)
            - `drug_col`  : drug name / pert_iname
            - `cell_col`  : LINCS cell line id
            - `dose_col`  : pert_dose (numeric, µM)
            - `time_col`  : pert_time (numeric, hours)
        Optionally also `is_gold_col` (0/1).

    cell_lines : None | str | iterable[str]
        - None: keep all cell lines (the "full sRGES" / SD6 variant).
        - str or iterable: restrict to those cell lines BEFORE the diff
          model is built — yields the per-cell-line / Figure 3 variant
          (e.g. cell_lines=["HEPG2"] for LIHC).

    cell_line_weights : None | Mapping[str, float]
        Per-cell-line `cor` values (e.g. correlation with target TCGA tumour
        type). If None, every cell line gets weight 1.0 — equivalent to
        skipping the cor*RGES weighting.

    apply_cor_weight : bool
        If False, skip the final (cor / max_cor) multiplication even when
        weights are provided.

    srges_diff_mode : {'binchen_r', 'corrected'}
        See module docstring. Default 'binchen_r' faithfully reproduces the
        upstream R code; 'corrected' fixes the integer-index mapping.

    return_full : bool
        If True (default), return all the columns the R pipeline emits:
        [drug_col, mean, n, median, sd, sRGES]. If False, return just
        [drug_col, sRGES].
    """
    if rges_df.empty:
        cols = [drug_col, "mean", "n", "median", "sd", "sRGES"] if return_full \
            else [drug_col, "sRGES"]
        return pd.DataFrame(columns=cols)

    df = rges_df.copy()
    
    # load LINCS metadata:
    lincs_metadata_df= pd.read_csv(LINCS_METADATA_PATH, usecols = ['id', cell_col,  dose_col, time_col], dtype='str')
    print(df.columns)
    # Add dose and time to drugs data
    df = df.merge(lincs_metadata_df,left_on="LINCS_id",  right_on="id", how="left")
    df = df.drop(columns=["id"])
    print(df.columns)


    # Optional is_gold filter.
    if filter_is_gold and is_gold_col is not None and is_gold_col in df.columns:
        df = df[df[is_gold_col].astype(str).isin(["1", "True", "TRUE"])]

    # Optional cell-line restriction.
    if cell_lines is not None:
        if isinstance(cell_lines, str):
            cell_lines = [cell_lines]
        df = df[df[cell_col].isin(list(cell_lines))]

    # Coerce numerics & drop unusable rows.
    for c in (score_col, dose_col, time_col):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=[score_col, dose_col, time_col, drug_col, cell_col])

    # ---- Step 1: build diff lookup
    diff = build_reference_diff(
        df,
        score_col=score_col, drug_col=drug_col, cell_col=cell_col,
        dose_col=dose_col, time_col=time_col,
        ref_dose=ref_dose, ref_time=ref_time,
        low_dose_threshold=low_dose_threshold,
        valid_times=valid_times,
    )
    diff_lookup = diff.to_dict()

    # ---- Step 2: per-instance adjustment
    # Only signatures with pert_dose > 0 and pert_time in valid_times survive,
    # matching the `pred` filter chain in sRGES_all_cmpds.R.
    pred = df[(df[dose_col] > 0) & (df[time_col].isin(valid_times))].copy()

    if cell_line_weights is None:
        pred["__cor__"] = 1.0
    else:
        pred["__cor__"] = pred[cell_col].map(cell_line_weights)
        # Bin Chen merges by cell_id and drops unmatched cell lines.
        pred = pred.dropna(subset=["__cor__"])

    if pred.empty:
        logger.warning(
            "compute_sRGES: no surviving rows after filtering "
            "(dose>0, time in %s, cor lookup). Returning empty.", valid_times,
        )
        cols = [drug_col, "mean", "n", "median", "sd", "sRGES"] if return_full \
            else [drug_col, "sRGES"]
        return pd.DataFrame(columns=cols)

    max_cor = pred["__cor__"].max()
    if max_cor == 0 or pd.isna(max_cor):
        max_cor = 1.0

    # Vectorize the row-wise function.
    rges_vals = pred[score_col].to_numpy()
    cor_vals = pred["__cor__"].to_numpy()
    dose_vals = pred[dose_col].to_numpy()
    time_vals = pred[time_col].to_numpy()

    adjusted = np.array([
        _get_sRGES_one(
            rges_vals[i], cor_vals[i], dose_vals[i], time_vals[i],
            diff_lookup, max_cor,
            low_dose_threshold=low_dose_threshold,
            short_time_threshold=short_time_threshold,
            apply_cor_weight=apply_cor_weight,
            srges_diff_mode=srges_diff_mode,
        )
        for i in range(len(pred))
    ])
    pred["__sRGES_instance__"] = adjusted

    # ---- Step 3: per-drug aggregation (sRGES = mean of adjusted instances)
    grouped = pred.groupby(drug_col, as_index=False).agg(
        mean=("__sRGES_instance__", "mean"),
        n=("__sRGES_instance__", "count"),
        median=("__sRGES_instance__", "median"),
        sd=("__sRGES_instance__", "std"),
    )
    grouped["sRGES"] = grouped["mean"]

    grouped = grouped.sort_values("sRGES").reset_index(drop=True)

    if return_full:
        return grouped[[drug_col, "mean", "n", "median", "sd", "sRGES"]]
    return grouped[[drug_col, "sRGES"]]


# ---------------------------------------------------------------------------
# Convenience: build cell-line weight map from BinChen-style files
# ---------------------------------------------------------------------------

def load_cell_line_weights_binchen(
    cell_line_cancer_csv: str,
    ccle_lincs_csv: str,
    *,
    primary_name_col: str = "Cell.line.primary.name",
    cor_col: str = "cor",
    ccle_name_col: str = "ccle_cell_line_name",
    lincs_id_col: str = "lincs_cell_id",
) -> dict:
    """
    Reproduce the merge step in sRGES_all_cmpds.R that maps the per-cell-line
    `cor` (TCGA correlation) to LINCS cell_id::

        cell_line_cancer <- merge(cell_line_cancer, ccle_lincs,
            by.x="Cell.line.primary.name", by.y="ccle_cell_line_name")

    Returns dict: lincs_cell_id -> cor.
    """
    cl = pd.read_csv(cell_line_cancer_csv)
    ccle = pd.read_csv(ccle_lincs_csv)
    merged = cl.merge(ccle, left_on=primary_name_col, right_on=ccle_name_col)
    return dict(zip(merged[lincs_id_col].astype(str), merged[cor_col].astype(float)))
