#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute correlations and classification metrics between a drug-ranking score
(e.g. RGES / sRGES / NetCoS connectivity score) and ChEMBL IC50 values for
MCF7, HepG2, and HT29 — the Bin Chen 2017 validation set.

Outputs
-------
- Spearman correlation between the drug-ranking score and the median IC50
  (linear and log10 scales).
- Precision / recall at user-defined CS and IC50 thresholds.
- Optional comparison against BinChen2017 SD5 sRGES values (same disease).
- Optional Bin Chen Fig.3-style scatter plot saved under `IMG_DIR`.
- Optional append of a one-row summary to `chembl_val_log_filename`.

Inputs (resolved via `src/conf.py` — no CLI)
--------------------------------------------
- CS rankings TSV produced by the NetCoS connectivity-score step, picked
  out of `logs/cs_runs.tsv` by `loader.pick_canonical_cs_run` from the
  upstream hyperparameters (CS_METHOD, cs_mith, cs_on_LM, CS_ON_PATHWAYS,
  landmark_disease, landmark_drug, cell_line).
- `BinChen2017/SD8.xlsx` — raw ChEMBL IC50 per disease (loaded by
  `loader.load_IC50`). Carries `pert_iname`, `pert_id`, `standard_value`,
  `standard_value_median`, `standard_inchi`, etc.
- `BinChen2017/SD5.xlsx` — BinChen sRGES + IC50 reference table per
  disease (used only when `compute_sd5=True`).

Merge semantics
---------------
Bin Chen 2017 joins drug-side tables on InChI (the canonical chemical
identifier), so this script follows the same convention via
`_merge_cs_ic50_inchi_first`:

  Pass 1 (primary): inner merge of CS rankings vs ChEMBL IC50 on
                    `standard_inchi` for rows whose InChI is non-empty on
                    both sides.
  Pass 2 (fallback): for rows not matched in pass 1, an inner merge on
                    `pert_iname`. Keeps coverage when InChI is missing on
                    either side; degenerates to the legacy single-pass
                    `pert_iname` merge when neither side carries InChI.

The SD5 overlap reporting (in the `compute_sd5` block) also prefers InChI
and falls back to `pert_iname` when SD5 doesn't expose `standard_inchi`.

Where `standard_inchi` comes from on each side
----------------------------------------------
- IC50 side: read directly from SD8 by `loader.load_IC50` when the column
  exists in the sheet.
- CS side: BinChen-replication shortcut: `LINCS_id` -> `pert_id` via the
  LINCS metadata file, then `pert_id` -> `standard_inchi` via SD8 itself.
  CS rows whose `pert_id` can't be resolved against SD8 are dropped. This
  path matches the strict BinChen2017 replication regime; once a richer
  external LINCS -> InChI mapping is available it should replace the SD8
  shortcut so we stop discarding CS rows that aren't in SD8.

Configuration
-------------
All runtime parameters (which CS run to load, which IC50 file, which cell
line, collapse methods, thresholds) come from `src/conf.py`. Swap
experiments by copying a `configs/<experiment>/conf.py` over `src/conf.py`
(see `configs/README.txt`). There is no CLI.

References
----------
- Bin Chen et al., Nature Communications (2017): "Reversal of cancer gene
  expression correlates with drug efficacy and reveals therapeutic
  targets." Validation correlates reversal potency with IC50 in MCF7
  (BRCA), HepG2 (LIHC), HT29 (COAD) via Spearman correlation against
  median IC50.
"""

import os
import sys
import warnings
HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))  # go up 3 levels
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

# Silence openpyxl's noisy "Unknown extension is not supported" warning that
# fires whenever pd.read_excel opens BinChen2017/{SD8,SD5}.xlsx. The files
# carry an Excel extension (custom XML namespace, conditional formatting, ...)
# that openpyxl's reader does not recognize; the actual data still loads
# correctly. Targeted filter — leaves all other openpyxl warnings intact.
warnings.filterwarnings(
    "ignore",
    message="Unknown extension is not supported",
    category=UserWarning,
    module="openpyxl",
)

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from conf import DISEASE, CS_OUT, DATA_DIR, CS_DIR,IMG_DIR,\
    cell_line, diseases_of, LOGS_DIR,\
        disease_run_name, cell_line_run_name,\
    cs_on_LM, cs_mith, selected_cs_run_id, \
    cs_log_filename, lincs_metadata_path, chembl_val_log_filename, ic50_file,\
    LINCS_METADATA_PATH, IC50_ONLY,CS_TH, IC50_EFF_TH, CS_DRUG_COLLAPSE_METHOD,\
        IC50_DRUG_COLLAPSE_METHOD, IC_50_binchen_SD5,\
    CS_METHOD, CS_ON_PATHWAYS

from logger import append_run_metadata

from loader import load_drug_rankings, load_IC50, pick_canonical_cs_run

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr, linregress, ttest_ind


def add_pert_iname_to_cs(cs_df, lincs_metadata_path, cs_id_col="LINCS_id",
                         metadata_id_col="id", metadata_name_col="pert_iname",
                         metadata_inchi_col="standard_inchi"):
    """Add `pert_iname` (and InChI when available) to a CS dataframe by
    left-joining LINCS metadata on LINCS_id.

    If `metadata_inchi_col` is present in the metadata file, it is also
    propagated so downstream code can merge CS vs IC50 on InChI (Bin Chen
    2017 convention). If the column is missing, only `pert_iname` is
    propagated; this matches the current BinChen-replication regime where
    InChI on the CS side is sourced from SD8 via
    `loader.add_pert_id_to_cs` + a downstream pert_id -> SD8 lookup rather
    than from `lincs_sig_info_new.csv`.
    """
    header_cols = pd.read_csv(lincs_metadata_path, nrows=0, dtype="str").columns
    use_cols = [metadata_id_col, metadata_name_col]
    if metadata_inchi_col in header_cols:
        use_cols.append(metadata_inchi_col)
    meta = pd.read_csv(lincs_metadata_path, usecols=use_cols, dtype="str")\
        .drop_duplicates(subset=[metadata_id_col])
    return cs_df.merge(meta, left_on=cs_id_col, right_on=metadata_id_col, how="left")\
        .drop(columns=[metadata_id_col])


def collapse_profiles_to_drug(in_df, score_col="connectivity_score", drug_col="pert_iname", how="median", keep_cols=None,
                              srges_kwargs=None):
    """Collapse multiple profiles per drug into one row per drug_col, preserving optional extra columns with 'first'.

    Parameters
    ----------
    how : str, one of {'median', 'mean', 'best', 'srges'}
        - 'median' / 'mean' / 'best' use a single-column groupby aggregation
          ('best' == min, since lower CS / lower RGES means stronger reversal).
        - 'srges' delegates to src/sRGES.py::compute_sRGES, which reproduces
          Bin Chen 2017's summarized RGES (with per-instance dose/time
          correction and optional cell-line cor weighting). This requires the
          input dataframe to also carry the per-signature metadata columns
          named in `srges_kwargs` (cell_col / dose_col / time_col,
          default 'cell_id' / 'pert_dose' / 'pert_time'). NetCoS CS values
          are fed directly 

    srges_kwargs : dict | None
        Extra keyword arguments forwarded to compute_sRGES when how='srges'.
        Common keys: cell_col, dose_col, time_col, cell_lines (str | list to
        restrict to specific LINCS cell lines — yields the Fig.3 per-cell-line
        variant; None yields the SD6 full-sRGES variant), cell_line_weights
        (dict mapping cell_id -> cor for the cor*RGES weighting; None to skip
        it), filter_is_gold, srges_diff_mode.
    """
    df = in_df.copy()
    if drug_col not in df.columns:
        raise ValueError(f"drug_col '{drug_col}' not found in dataframe")
    if score_col not in df.columns:
        raise ValueError(f"score_col '{score_col}' not found in dataframe")

    keep_cols = [] if keep_cols is None else keep_cols
    missing_keep_cols = [col for col in keep_cols if col not in df.columns]
    if missing_keep_cols:
        raise ValueError(f"keep_cols not found in dataframe: {missing_keep_cols}")

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[drug_col, score_col]).copy()

    if how == "srges":
        # Lazy import to keep sRGES.py optional for callers that never use it.
        from sRGES import compute_sRGES
        srges_kwargs = dict(srges_kwargs or {})
        out = compute_sRGES(df, score_col=score_col, drug_col=drug_col,
            return_full=False, **srges_kwargs)
        # Match the return shape of the median/mean/best branch: drug_col +
        # score_col (with sRGES value written into score_col), plus keep_cols.
        out = out.rename(columns={"sRGES": score_col})
        if keep_cols:
            first = df.groupby(drug_col, as_index=False)[keep_cols].first()
            out = out.merge(first, on=drug_col, how="left")
        return out

    if how == "median":
        score_agg = "median"
    elif how == "mean":
        score_agg = "mean"
    elif how == "best":
        score_agg = "min"
    else:
        raise ValueError("how must be one of: 'median', 'mean', 'best', 'srges'")

    agg_dict = {score_col: score_agg, **{col: "first" for col in keep_cols}}
    return df.groupby(drug_col, as_index=False).agg(agg_dict)

def _coerce_cs_on_pathways(v):
    """Normalize CS_ON_PATHWAYS values from cs_runs.tsv (bool / int / str) to bool."""
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return bool(int(v))
    return str(v).strip().lower() in ('true', '1', '1.0')


def resolve_cs_run_id(cs_runs_tsv, disease_run_id, drug_run_id, cs_on_LM, mith,
                      CS_METHOD=None, CS_ON_PATHWAYS=None, selected_cs_run_id=None):
    """
    Retrieve the cs_run_id from cs log file, or return selected_cs_run_id, if provided.

    CS_METHOD and CS_ON_PATHWAYS are optional filters (default None = no filter).
    They are needed for grid-search use cases where cs_runs.tsv contains multiple
    rows sharing (disease_run_id, drug_run_id, cs_on_LM, mith) but differing in
    method / pathways. Old callers that don't pass them keep their previous
    "most-recent matching row" semantics.
    """

    runs = pd.read_csv(cs_runs_tsv, sep="\t")

    # normalize booleans stored as 0/1 or True/False
    runs["mith"] = runs["mith"].astype(int)
    runs["cs_on_LM"] = runs["cs_on_LM"].astype(int)

    if selected_cs_run_id is not None:
        hit = runs[runs["cs_run_id"] == selected_cs_run_id]
        if len(hit) == 0:
            raise ValueError(
                f"selected_cs_run_id '{selected_cs_run_id}' not found in {cs_runs_tsv}"
            )
        return selected_cs_run_id

    mask = (
        (runs["disease_run_id"] == disease_run_id) &
        (runs["drug_run_id"] == drug_run_id) &
        (runs["cs_on_LM"] == int(cs_on_LM)) &
        (runs["mith"] == int(mith))
    )
    if CS_METHOD is not None:
        mask &= (runs["CS_METHOD"] == CS_METHOD)
    if CS_ON_PATHWAYS is not None:
        mask &= runs["CS_ON_PATHWAYS"].apply(_coerce_cs_on_pathways) == bool(CS_ON_PATHWAYS)

    hit = runs[mask].copy()

    if len(hit) == 0:
        raise ValueError(
            "No matching CS run found in cs_runs.tsv for:\n"
            f"disease_run_id={disease_run_id}, "
            f"drug_run_id={drug_run_id}, "
            f"cs_on_LM={int(cs_on_LM)}, "
            f"mith={int(mith)}, "
            f"CS_METHOD={CS_METHOD}, "
            f"CS_ON_PATHWAYS={CS_ON_PATHWAYS}"
        )

    if len(hit) > 1:
        # cs_runs.tsv mixes dd/mm/yyyy HH:MM (older rows) with ISO-8601
        # (newer rows). dayfirst=True parses the dd/mm form correctly and
        # is ignored for the unambiguous ISO form.
        hit["timestamp"] = pd.to_datetime(hit["timestamp"], dayfirst=True)
        hit = hit.sort_values("timestamp", ascending=False)

    return hit.iloc[0]["cs_run_id"]

# -----------------------------------------------------------------------------
# Merge helpers (CS rankings <-> ChEMBL IC50)
# -----------------------------------------------------------------------------

def _coalesce_suffixed(merged_df, base_col, suffixes=("_dr", "_ic50")):
    """If both `base_col + suffixes[0]` and `base_col + suffixes[1]` exist
    in merged_df, collapse them into a single `base_col` (first non-null
    wins) and drop the suffixed pair. Used to normalize column schemas
    across the InChI-keyed and pert_iname-keyed merge passes so they can be
    concatenated cleanly.
    """
    a = base_col + suffixes[0]
    b = base_col + suffixes[1]
    if a in merged_df.columns and b in merged_df.columns:
        merged_df = merged_df.copy()
        merged_df[base_col] = merged_df[a].combine_first(merged_df[b])
        merged_df = merged_df.drop(columns=[a, b])
    return merged_df


def _merge_cs_ic50_inchi_first(dr, ic50, *, cs_drug_colname='pert_iname',
                               ic50_drug_colname='pert_iname',
                               inchi_col='standard_inchi', verbose=False):
    """Two-pass inner merge of CS rankings against IC50 measurements.

    Bin Chen 2017 uses InChI as the canonical chemical identifier when
    joining drug-side tables, so this merge prefers InChI over pert_iname.

    Pass 1 (primary): inner merge on `inchi_col` for rows whose InChI is
        non-empty on both sides.
    Pass 2 (fallback): inner merge on `pert_iname` for rows not matched in
        pass 1. Covers compounds with missing InChI on either side (LINCS
        metadata not yet augmented, ChEMBL entries with empty standard_inchi,
        etc.) so we don't lose coverage during the transition.

    The two pass results are normalized to the same column schema before
    being concatenated. The returned dataframe carries single un-suffixed
    `pert_iname` and `standard_inchi` columns and behaves like the original
    single-pass merge result for all downstream consumers.
    """
    use_inchi = (
        inchi_col in dr.columns and inchi_col in ic50.columns
    )

    if not use_inchi:
        # Fall back to the legacy single-pass merge.
        merged = pd.merge(
            dr, ic50,
            left_on=cs_drug_colname, right_on=ic50_drug_colname,
            how='inner', suffixes=('_dr', '_ic50'),
        )
        if verbose:
            print(f"[merge] InChI column '{inchi_col}' not available on both "
                  f"sides; falling back to single-pass merge on "
                  f"{cs_drug_colname}.")
        return merged

    # --- Pass 1: InChI-keyed merge ---
    dr_inchi_mask = dr[inchi_col].notna() & (dr[inchi_col].astype(str).str.strip() != '')
    ic50_inchi_mask = ic50[inchi_col].notna() & (ic50[inchi_col].astype(str).str.strip() != '')
    dr_with_inchi = dr[dr_inchi_mask].copy()
    ic50_with_inchi = ic50[ic50_inchi_mask].copy()

    if not dr_with_inchi.empty and not ic50_with_inchi.empty:
        merged_inchi = pd.merge(
            dr_with_inchi, ic50_with_inchi,
            on=inchi_col, how='inner',
            suffixes=('_dr', '_ic50'),
        )
        # After merging on InChI, both sides retain their pert_iname
        # column (suffixed). Coalesce to a single pert_iname.
        merged_inchi = _coalesce_suffixed(merged_inchi, cs_drug_colname)
    else:
        merged_inchi = pd.DataFrame()

    # --- Pass 2: pert_iname fallback ---
    if cs_drug_colname in merged_inchi.columns and not merged_inchi.empty:
        matched_inames = set(merged_inchi[cs_drug_colname].dropna().astype(str))
    else:
        matched_inames = set()

    dr_left = dr[~dr[cs_drug_colname].astype(str).isin(matched_inames)].copy()
    ic50_left = ic50[~ic50[ic50_drug_colname].astype(str).isin(matched_inames)].copy()

    if not dr_left.empty and not ic50_left.empty:
        merged_name = pd.merge(
            dr_left, ic50_left,
            left_on=cs_drug_colname, right_on=ic50_drug_colname,
            how='inner', suffixes=('_dr', '_ic50'),
        )
        # pert_iname is single (left_on==right_on collapses it); InChI gets
        # suffixed if both sides have it -> coalesce.
        merged_name = _coalesce_suffixed(merged_name, inchi_col)
    else:
        merged_name = pd.DataFrame()

    if verbose:
        print(f"[merge] inchi-matched rows: {len(merged_inchi)}; "
              f"pert_iname-fallback rows: {len(merged_name)}")

    merged = pd.concat([merged_inchi, merged_name], ignore_index=True, sort=False)
    return merged


# classification and plotting

def classify_ic50_vs_cs(df, score_col="connectivity_score", ic50_col="standard_value_median", cs_threshold = -1.5, ic50_threshold = 10.0,):
    """
    Classify compounds into TP / FP / TN / FN.

    Positive class = effective drug.
    - predicted positive: CS <= cs_threshold
      (more negative CS = stronger predicted reversal)
    - actual positive: IC50 <= ic50_threshold
      (lower IC50 = stronger efficacy)
    """
    out = df.copy()
    out[score_col] = pd.to_numeric(out[score_col], errors="coerce")
    out[ic50_col] = pd.to_numeric(out[ic50_col], errors="coerce")
    out = out.dropna(subset=[score_col, ic50_col]).copy()

    out["predicted_positive"] = out[score_col] <= float(cs_threshold)
    out["actual_positive"] = out[ic50_col] <= float(ic50_threshold)

    out["quadrant"] = np.select(
        [
            out["predicted_positive"] & out["actual_positive"],
            out["predicted_positive"] & ~out["actual_positive"],
            ~out["predicted_positive"] & out["actual_positive"],
            ~out["predicted_positive"] & ~out["actual_positive"],
        ],
        ["TP", "FP", "FN", "TN"],
        default="NA",
    )

    tp = int((out["quadrant"] == "TP").sum())
    fp = int((out["quadrant"] == "FP").sum())
    fn = int((out["quadrant"] == "FN").sum())
    tn = int((out["quadrant"] == "TN").sum())

    precision = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    recall = tp / (tp + fn) if (tp + fn) > 0 else np.nan

    metrics = {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": precision,
        "recall": recall,
        "cs_threshold": float(cs_threshold),
        "ic50_threshold": float(ic50_threshold),
    }
    return out, metrics

def select_top_residual_annotations(df, x, y, lr, drug_label_col="drug", annotate_top_n=5, min_dx_frac=0.03, min_dy_frac=0.03,):
    """
    Select rows to annotate based on largest absolute residuals from the
    regression line, while avoiding repeated labels and points that are too
    close to each other.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame corresponding exactly to x and y.
    x : np.ndarray
        X coordinates used in the scatter plot.
    y : np.ndarray
        Y coordinates used in the scatter plot.
    lr : scipy.stats._stats_py.LinregressResult
        Result of scipy.stats.linregress(x, y).
    drug_label_col : str
        Column containing the label to annotate.
    annotate_top_n : int
        Maximum number of annotations to return.
    min_dx_frac : float
        Minimum allowed x-distance between annotated points, as a fraction
        of total x range.
    min_dy_frac : float
        Minimum allowed y-distance between annotated points, as a fraction
        of total y range.

    Returns
    -------
    selected_df : pd.DataFrame
        Subset of rows selected for annotation, with helper columns:
        '_x', '_y', '_resid', '_label'.
    """
    if annotate_top_n is None or annotate_top_n <= 0:
        return pd.DataFrame(columns=list(df.columns) + ["_x", "_y", "_resid", "_label"])

    if drug_label_col not in df.columns:
        return pd.DataFrame(columns=list(df.columns) + ["_x", "_y", "_resid", "_label"])

    if len(df) != len(x) or len(df) != len(y):
        raise ValueError("df, x, and y must have the same length.")

    work_df = df.copy()
    work_df["_x"] = np.asarray(x, dtype=float)
    work_df["_y"] = np.asarray(y, dtype=float)
    work_df["_label"] = work_df[drug_label_col].astype(str)

    pred = lr.intercept + lr.slope * work_df["_x"].to_numpy()
    work_df["_resid"] = np.abs(work_df["_y"].to_numpy() - pred)

    # Sort from largest residual to smallest
    work_df = work_df.sort_values("_resid", ascending=False).reset_index(drop=True)

    x_range = work_df["_x"].max() - work_df["_x"].min()
    y_range = work_df["_y"].max() - work_df["_y"].min()

    # Fallback values for degenerate ranges
    min_dx = min_dx_frac * x_range if x_range > 0 else 0.01
    min_dy = min_dy_frac * y_range if y_range > 0 else 0.01

    used_labels = set()
    used_positions = []
    selected_idx = []

    for idx, row in work_df.iterrows():
        label = row["_label"]
        px = row["_x"]
        py = row["_y"]

        # Skip repeated labels
        if label in used_labels:
            continue

        # Skip rows too close to an already selected annotation point
        too_close = any(
            abs(px - qx) < min_dx and abs(py - qy) < min_dy
            for qx, qy in used_positions
        )
        if too_close:
            continue

        selected_idx.append(idx)
        used_labels.add(label)
        used_positions.append((px, py))

        if len(selected_idx) >= annotate_top_n:
            break

    return work_df.loc[selected_idx].copy()

def annotate_selected_points(ax, selected_df, label_col="_label", x_col="_x", y_col="_y", fontsize=9, xytext=(4, 4), add_arrows=False,):
    """
    Annotate selected points on an axis.
    """
    for _, row in selected_df.iterrows():
        kwargs = {
            "fontsize": fontsize,
            "xytext": xytext,
            "textcoords": "offset points",
        }
        if add_arrows:
            kwargs["arrowprops"] = dict(arrowstyle="-", lw=0.5)

        ax.annotate(
            str(row[label_col]),
            (row[x_col], row[y_col]),
            **kwargs,
        )

def plot_binchen_fig3_style(
    classified_df,
    metrics,
    disease,
    cell_line_name,
    score_col="connectivity_score",
    ic50_col="standard_value",
    log_ic50_col="log10_ic50",
    drug_label_col="pert_iname",
    quadrant_col="quadrant",
    annotate_top_n=5,
    output_file=None,
    cs_threshold=-1.5,
    ic50_threshold=10.0,
):
    """
    Reproduce the layout of Bin Chen Fig. 3:
        - scatter of score vs log10(IC50)
        - fitted regression line
        - annotations for compounds with largest residuals
        - boxplot of score in effective vs ineffective compounds
            and adds a TP/FP/FN/TN threshold cross.

    Returns
    -------
    fig, axes, stats_dict
    """
    df = classified_df.copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df[ic50_col] = pd.to_numeric(df[ic50_col], errors="coerce")
    df = df.dropna(subset=[score_col, ic50_col]).copy()
    df[log_ic50_col] = np.log10(df[ic50_col])

    if log_ic50_col not in df.columns:
        df[log_ic50_col] = np.log10(df[ic50_col])

    x = df[score_col].to_numpy(dtype=float)
    y = df[log_ic50_col].to_numpy(dtype=float)

    lr = linregress(x, y)
    rho, rho_p = spearmanr(x, y)

    effective_scores = classified_df.loc[
        classified_df[ic50_col] <= ic50_threshold, score_col
    ].astype(float).to_numpy()
    ineffective_scores = classified_df.loc[
        classified_df[ic50_col] > ic50_threshold, score_col
    ].astype(float).to_numpy()

    t_res = ttest_ind(effective_scores, ineffective_scores, equal_var=False)

    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(12, 5),
        gridspec_kw={"width_ratios": [3.2, 1.2]},
    )

    # --- left panel: scatter ---
    ax = axes[0]

    # plot regression line
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 200)
    yy = lr.intercept + lr.slope * xx
    ax.plot(xx, yy, color='darkgrey')
    
    # scatter plot    

    ax.scatter(x, y, s=8, facecolors= 'none', \
               edgecolors=np.where(df[ic50_col] <= ic50_threshold, "cornflowerblue", "lightcoral"), lw=0.5)


    ax.set_xlabel("Connectivity score")
    ax.set_ylabel("Log10 (IC50)")
    ax.set_title(f"{disease}, {cell_line_name}")

    # threshold cross
    y_thr = np.log10(ic50_threshold)
    ax.axvline(cs_threshold, linestyle="--", linewidth=1, color='red', alpha=0.2)
    ax.axhline(y_thr, linestyle="--", linewidth=1, color='red', alpha=0.2)
    
    
    
    # add legend
    txt = (
        f"r={lr.rvalue:.2f}, P={lr.pvalue:.2e}\n"
        f"rho={rho:.2f}, P={rho_p:.2e}\n"
        f"precision={metrics['precision']:.2f}, recall={metrics['recall']:.2f}"
    )
    ax.text(0.625,0.97, txt, transform=ax.transAxes, va="top",ha="left",  bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),  )
    
    # annotate true positive drugs
    tp_mask = (df[score_col] <= cs_threshold) & (df[ic50_col] <= ic50_threshold)

    for _, row in df.loc[tp_mask].iterrows():
        ax.annotate(str(row[drug_label_col]),(row[score_col], row[log_ic50_col]),fontsize=8, xytext=(4, 4),textcoords="offset points",)
    
    # annotate largest residuals, as in Bin Chen figure
    if annotate_top_n and drug_label_col in df.columns:
        selected_annotations = select_top_residual_annotations(df=df, x=x, y=y,\
                                lr=lr,drug_label_col=drug_label_col,annotate_top_n=annotate_top_n, min_dx_frac=0.03,min_dy_frac=0.03,)

        annotate_selected_points(ax=ax, selected_df=selected_annotations,fontsize=9,xytext=(4, 4),add_arrows=False,)

    # quadrant labels
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    
    # TP quadrant highlighting
    tp_rect = Rectangle(
        (x_min, y_min),                     # bottom-left corner
        cs_threshold - x_min,               # width
        y_thr - y_min,                      # height
        facecolor="lightblue",
        alpha=0.2,
        zorder=0
    )
    ax.add_patch(tp_rect)

    def _mid(a, b):
        return a + (b - a) / 2.0

    quadrant_specs = [
        (_mid(x_min, cs_threshold), _mid(y_min, y_thr), f"TP\nn={metrics['tp']}"),
        (_mid(x_min, cs_threshold), _mid(y_thr, y_max), f"FP\nn={metrics['fp']}"),
        (_mid(cs_threshold, x_max), _mid(y_min, y_thr), f"FN\nn={metrics['fn']}"),
        (_mid(cs_threshold, x_max), _mid(y_thr, y_max), f"TN\nn={metrics['tn']}"),
    ]
    for qx, qy, label in quadrant_specs:
        ax.text(qx, qy, label, ha="center", va="center", alpha=0.75)

    # --- right panel: boxplot ---
    ax2 = axes[1]
    box_data = [effective_scores, ineffective_scores]
    bp=ax2.boxplot(box_data, tick_labels=["Effective", "Ineffective"],  medianprops=dict(color="black", linewidth=2), patch_artist=True)
    [b.set_facecolor(c) for b, c in zip(bp["boxes"], ["cornflowerblue", "lightcoral"])]
    ax2.set_ylabel("Connectivity score")
    ax2.set_title(f"P={t_res.pvalue:.2e}")

    fig.tight_layout()

    if output_file is not None:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    stats_dict = {
        "n_compounds": int(classified_df.shape[0]),
        "n_effective": int((classified_df[ic50_col] <= ic50_threshold).sum()),
        "n_ineffective": int((classified_df[ic50_col] > ic50_threshold).sum()),
        "pearson_r": float(lr.rvalue),
        "pearson_p": float(lr.pvalue),
        "spearman_rho": float(rho),
        "spearman_p": float(rho_p),
        "ttest_p": float(t_res.pvalue),
        "tp": metrics["tp"],
        "fp": metrics["fp"],
        "fn": metrics["fn"],
        "tn": metrics["tn"],
        "precision": float(metrics["precision"]) if pd.notna(metrics["precision"]) else np.nan,
        "recall": float(metrics["recall"]) if pd.notna(metrics["recall"]) else np.nan,
        "cs_threshold": float(cs_threshold),
        "ic50_threshold": float(ic50_threshold),
    }

    return fig, axes, stats_dict


# -----------------------------------------------------------------------------
# Reusable per-conf routine
# -----------------------------------------------------------------------------

def run_correlation_for_conf(
    *,
    disease,
    cell_line,
    landmark_disease,
    landmark_drug,
    CS_METHOD,
    cs_mith,
    cs_on_LM,
    CS_ON_PATHWAYS,
    cs_drug_collapse_method,
    ic50_drug_collapse_method,
    ic50_only,
    cell_line_IC50,
    cs_threshold=-1.5,
    ic50_threshold=10.0,
    selected_cs_run_id=None,
    compute_sd5=False,
    plot=False,
    log_to_tsv=False,
    annotate_top_n=5,
    verbose=True,
):
    """
    Run the full ChEMBL IC50 vs CS validation routine for one set of
    hyperparameters. Equivalent to the legacy __main__ block, but parameterized.

    Merge semantics
    ---------------
    CS rankings and ChEMBL IC50 are joined via `_merge_cs_ic50_inchi_first`:
    primary inner merge on `standard_inchi`, fallback inner merge on
    `pert_iname` for rows with missing InChI on either side. The SD5
    overlap (when `compute_sd5=True`) is computed on the same key with the
    same fallback. See the module-level docstring for the full description.

    Parameters
    ----------
    cell_line, landmark_disease, landmark_drug, CS_METHOD, cs_mith, cs_on_LM,
    CS_ON_PATHWAYS
        Pre-MITHrIL / pre-CS hyperparameters that identify a unique row in
        cs_runs.tsv (and therefore a unique CS output file).

    cs_drug_collapse_method, ic50_drug_collapse_method, ic50_only,
    cell_line_IC50, cs_threshold, ic50_threshold
        Validation hyperparameters applied on top of the resolved CS file.
        cell_line_IC50 may be a cell-line string or None ("use all cell
        lines in the IC50 file").

    selected_cs_run_id : optional override. If provided, skips resolution.
    compute_sd5 : if True, also compute spearman / precision / recall against
        the BinChen 2017 SD5 sRGES table for the same disease. The
        InChI-vs-pert_iname overlap key actually used is logged in the
        returned summary dict under `overlap_key_BC_SD5`.
    plot : if True, save a Bin Chen Fig.3-style plot to IMG_DIR.
    log_to_tsv : if True, append the summary dict as a row in
        chembl_val_log_filename.

    Returns
    -------
    classified_df : pd.DataFrame
    pr_metrics : dict
    summary : dict
    """
    cs_drug_colname = 'pert_iname'
    ic50_drug_colname = 'pert_iname'
    ic50_score_col = 'standard_value'

    # ---- resolve identifiers from upstream params ----
    disease = diseases_of[cell_line]
    disease_run_id = disease + ('_LM' if landmark_disease else '')
    drug_run_id = cell_line + ('_LM' if landmark_drug else '')
    cs_out_dir = CS_DIR / 'output' / disease_run_id

    if verbose:
        print('running correlations between chembl IC50 and drug rankings for disease:', disease)
        print(disease, cs_out_dir, cell_line)

    # cs_run_id = resolve_cs_run_id(
    #     cs_runs_tsv=cs_log_filename,
    #     disease_run_id=disease_run_id,
    #     drug_run_id=drug_run_id,
    #     cs_on_LM=cs_on_LM,
    #     mith=cs_mith,
    #     CS_METHOD=CS_METHOD,
    #     CS_ON_PATHWAYS=CS_ON_PATHWAYS,
    #     selected_cs_run_id=selected_cs_run_id,
    # )
    # cs_drug_file = f"{cs_run_id}.tsv"
    
    cs_drug_file = pick_canonical_cs_run(
        disease=disease,
        landmark_disease_pref=landmark_disease,
        cs_log_path=None,
        cs_dir=None,
        cs_method=CS_METHOD,
        mith=cs_mith,
        cs_on_LM=cs_on_LM,
        cs_on_pathways=CS_ON_PATHWAYS,
        require_file_on_disk=True,
        verbose=True,
    )
    
    cs_run_id = cs_drug_file.split('.')[0]
    
    if verbose:
        print("Using CS file:", cs_out_dir / cs_drug_file)

    # Bin Chen 2017 joins drug-side tables on InChI. The CS side picks it
    # up via add_pert_iname_to_cs (when the LINCS metadata file carries
    # standard_inchi), and the IC50 side picks it up via loader.load_IC50.
    inchi_col = 'standard_inchi'

    # ---- load CS rankings ----
    dr_uncollapsed = load_drug_rankings(cs_out_dir, filename=cs_drug_file)
    dr_uncollapsed = add_pert_iname_to_cs(dr_uncollapsed, LINCS_METADATA_PATH)
    if verbose:
        print(dr_uncollapsed.shape, 'drugs ')
    dr_keep = [inchi_col] if inchi_col in dr_uncollapsed.columns else []
    dr = collapse_profiles_to_drug(
        in_df=dr_uncollapsed,
        drug_col=cs_drug_colname,
        how=cs_drug_collapse_method,
        keep_cols=dr_keep,
    )
    if verbose:
        print(dr.shape, 'drugs ')

    # ---- load IC50 ----
    ic50_uncollapsed, ic50_log_data = load_IC50(
        ic50_file, disease,
        cell_line=cell_line_IC50,
        IC50_ONLY=ic50_only,
        IC_50_binchen_SD5=IC_50_binchen_SD5,
    )
    if verbose:
        print(ic50_uncollapsed.shape, 'drugs with IC50 value for cell line', cell_line_IC50)
    ic50_keep = [inchi_col] if inchi_col in ic50_uncollapsed.columns else []
    ic50 = collapse_profiles_to_drug(
        in_df=ic50_uncollapsed,
        score_col=ic50_score_col,
        drug_col=ic50_drug_colname,
        how=ic50_drug_collapse_method,
        keep_cols=ic50_keep,
    )
    if verbose:
        print(ic50.shape, 'drugs with IC50 value for cell line', cell_line_IC50)

    # ---- merge: InChI primary, pert_iname fallback ----
    merged = _merge_cs_ic50_inchi_first(
        dr, ic50,
        cs_drug_colname=cs_drug_colname,
        ic50_drug_colname=ic50_drug_colname,
        inchi_col=inchi_col,
        verbose=verbose,
    )
    merged = merged.dropna(subset=["connectivity_score", ic50_score_col]).copy()
    if verbose:
        print('merged data on common drugs:', merged.shape)
    merged["log10_ic50"] = np.log10(merged[ic50_score_col])

    # ---- spearman ----
    rho_linear, p_linear = spearmanr(merged["connectivity_score"], merged[ic50_score_col])
    rho_log, p_log = spearmanr(merged["connectivity_score"], merged["log10_ic50"])
    if verbose:
        print('linear IC50: rho=', np.round(rho_linear, 2), ' pval =', np.round(p_linear, 2))
        print('log IC50: rho=', np.round(rho_log, 2), ' pval =', np.round(p_log, 2))

    # ---- precision / recall ----
    classified_df, pr_metrics = classify_ic50_vs_cs(
        merged,
        score_col="connectivity_score",
        ic50_col=ic50_score_col,
        cs_threshold=cs_threshold,
        ic50_threshold=ic50_threshold,
    )
    if verbose:
        print(
            "classification:",
            f"TP={pr_metrics['tp']}, FP={pr_metrics['fp']}, "
            f"FN={pr_metrics['fn']}, TN={pr_metrics['tn']}, "
            f"precision={np.round(pr_metrics['precision'], 2)}, "
            f"recall={np.round(pr_metrics['recall'], 2)}",
        )

    # ---- optional SD5 (BinChen 2017 sRGES) overlap and metrics ----
    sd5_block = {}
    if compute_sd5:
        BC_merged = pd.read_excel(
            DATA_DIR / 'BinChen2017' / 'SD5.xlsx',
            sheet_name=disease,
        )
        # Overlap with the NetCoS merged result: join on InChI (BinChen
        # 2017 convention) when both sides carry `standard_inchi`. Falls
        # back to pert_iname if SD5 doesn't expose InChI, so this also
        # works on older SD5 sheets. The overlap is a sanity-check count;
        # the spearman below is computed on BC_merged alone (SD5's own
        # sRGES against SD5's IC50) and is independent of the join key.
        if inchi_col in BC_merged.columns and inchi_col in merged.columns:
            overlap_key = inchi_col
            cs_inchis = set(
                merged[inchi_col].dropna().astype(str).str.strip()
            ) - {''}
            sd5_inchis = set(
                BC_merged[inchi_col].dropna().astype(str).str.strip()
            ) - {''}
            overlap_BC = cs_inchis.intersection(sd5_inchis)
        else:
            overlap_key = cs_drug_colname
            overlap_BC = set(merged[cs_drug_colname]).intersection(
                set(BC_merged[cs_drug_colname])
            )
        if verbose:
            print(f"SD5 overlap on '{overlap_key}': {len(overlap_BC)} drugs")
            print(overlap_BC)
        BC_merged["log10_ic50"] = np.log10(BC_merged[ic50_score_col])
        rho_linear_bc, p_linear_bc = spearmanr(BC_merged["sRGES"], BC_merged[ic50_score_col])
        rho_log_bc, p_log_bc = spearmanr(BC_merged["sRGES"], BC_merged["log10_ic50"])
        if verbose:
            print('linear IC50 SD5: rho=', np.round(rho_linear_bc, 2), ' pval =', np.round(p_linear_bc, 2))
            print('log IC50 SD5: rho=', np.round(rho_log_bc, 2), ' pval =', np.round(p_log_bc, 2))
        sd5_block = {
            "spearman_r_SD5": np.round(rho_linear_bc, 2),
            "spearman_pval_SD5": np.round(p_linear_bc, 2),
            "spearman_r_log_SD5": np.round(rho_log_bc, 2),
            "spearman_pval_log_SD5": np.round(p_log_bc, 2),
            "n_overlap_BC_SD5": len(overlap_BC),
            "overlap_key_BC_SD5": overlap_key,
        }

    # ---- optional plot ----
    plot_file = None
    if plot:
        plot_file = IMG_DIR / (
            f"{disease}_{cs_run_id}_{cs_drug_collapse_method}_"
            f"{ic50_drug_collapse_method}_fig3_style.png"
        )
        plot_binchen_fig3_style(
            classified_df=classified_df,
            metrics=pr_metrics,
            disease=disease,
            cell_line_name=cell_line,
            score_col="connectivity_score",
            ic50_col=ic50_score_col,
            drug_label_col=cs_drug_colname,
            annotate_top_n=annotate_top_n,
            output_file=plot_file,
            cs_threshold=cs_threshold,
            ic50_threshold=ic50_threshold,
        )

    # ---- summary (one row in chembl_val_log_filename) ----
    summary = {
        # identity
        "IC50_validation_run_id": datetime.now().strftime("%d_%m_%Y_%H_%M_%S"),
        "datetime": datetime.now().isoformat(),

        # inputs
        "cs_file": str(cs_out_dir / cs_drug_file),
        "ic50_file": str(ic50_file),
        "lincs_metadata_file": str(lincs_metadata_path),
        "disease_run_id": disease_run_id,
        "drug_run_id": drug_run_id,
        "cs_run_id": cs_run_id,

        # pre-CS hyperparameters
        "cell_line": cell_line,
        "landmark_disease": landmark_disease,
        "landmark_drug": landmark_drug,
        "CS_METHOD": CS_METHOD,
        "cs_mith": cs_mith,
        "cs_on_LM": cs_on_LM,
        "CS_ON_PATHWAYS": CS_ON_PATHWAYS,

        # validation hyperparameters
        "cs_drug_colname": cs_drug_colname,
        "ic50_drug_colname": ic50_drug_colname,
        "LINCS_drug_collapse_method": cs_drug_collapse_method,
        "IC50_drug_collapse_method": ic50_drug_collapse_method,
        "ic50_only": ic50_only,
        "cell_line_IC50": cell_line_IC50,
        "cs_threshold": cs_threshold,
        "ic50_threshold": ic50_threshold,

        # sizes
        "n_cs_rows_before_drug_collapse": len(dr_uncollapsed),
        "n_cs_rows_after_drug_collapse": len(dr),
        "n_ic50_rows_before_drug_collapse": len(ic50_uncollapsed),
        "n_ic50_rows_after_drug_collapse": len(ic50),
        "n_merged_rows": len(merged),
        "n_unique_drugs": merged[cs_drug_colname].nunique(),

        # spearman
        "spearman_r": np.round(rho_linear, 2),
        "spearman_pval": np.round(p_linear, 2),
        "spearman_r_log": np.round(rho_log, 2),
        "spearman_pval_log": np.round(p_log, 2),

        # precision/recall
        "tp": pr_metrics["tp"],
        "fp": pr_metrics["fp"],
        "fn": pr_metrics["fn"],
        "tn": pr_metrics["tn"],
        "precision": (np.round(pr_metrics["precision"], 4)
                      if pd.notna(pr_metrics["precision"]) else np.nan),
        "recall": (np.round(pr_metrics["recall"], 4)
                   if pd.notna(pr_metrics["recall"]) else np.nan),

        # output
        "output_plot_file": str(plot_file) if plot_file is not None else None,
    }
    summary.update(sd5_block)
    summary.update(ic50_log_data)

    if log_to_tsv:
        append_run_metadata(chembl_val_log_filename, summary)

    return classified_df, pr_metrics, summary



if __name__ == "__main__":
    # Run the routine using conf.py defaults — preserves the original
    # single-disease interactive behavior.

    # The landmark flags are derived from the run-name suffixes that conf.py
    # already computed. Keeps a single source of truth.
    landmark_disease_flag = disease_run_name.endswith('_LM')
    landmark_drug_flag = cell_line_run_name.endswith('_LM')

    classified_df, pr_metrics, summary = run_correlation_for_conf(
        disease=DISEASE,
        cell_line=cell_line,
        landmark_disease=landmark_disease_flag,
        landmark_drug=landmark_drug_flag,
        CS_METHOD=CS_METHOD,
        cs_mith=cs_mith,
        cs_on_LM=cs_on_LM,
        CS_ON_PATHWAYS=CS_ON_PATHWAYS,
        cs_drug_collapse_method=CS_DRUG_COLLAPSE_METHOD,
        ic50_drug_collapse_method=IC50_DRUG_COLLAPSE_METHOD,
        ic50_only=IC50_ONLY,
        cell_line_IC50=cell_line,
        cs_threshold=-0.1,
        ic50_threshold=IC50_EFF_TH,
        selected_cs_run_id=selected_cs_run_id,
        compute_sd5=True,
        plot=True,
        log_to_tsv=True,
    )


