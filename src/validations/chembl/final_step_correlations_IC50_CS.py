#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute correlations and classification metrics between a drug-ranking score (e.g., RGES/sRGES/NetCos)
and ChEMBL IC50 values for MCF7, HepG2, and HT29 (Bin Chen 2017 validation).

-Spearman correlation
- Precision/recall of IC50 vs drug ranking score

Inputs
------
1) Drug ranking file (CSV/TSV): must contain at least two columns
   - "compound" (drug name or InChIKey)
   - "score"    (your ranking value; more negative/positive as per your method)
   Optional columns are ignored.

2) IC50 file(s) from Bin Chen 2017 Supplementary Data (CSV or XLSX).
   - If XLSX: the file may have one or multiple sheets; the script will read all sheets
     and look for columns that map to ["compound", "cell_line", "ic50"].
   - If CSV/TSV: same expectation.
   - If you have both Supplementary Data 3 and 8, you can pass multiple files; they will be concatenated.

Assumptions
-----------
- IC50 units are consistent within the file(s) (often µM in the Supplementary Data).
- Multiple IC50 measurements per (compound, cell_line) are reduced by median, as in the paper.
- Join key is "compound" (case-insensitive). If both datasets have an "inchi_key" column,
  set --join-key inchi_key to join on that instead.

Outputs
-------
- Prints a small summary table with N, Spearman rho, and p-value for each cell line.
- Optionally writes the merged tables per cell line to CSV for inspection (--out-dir).

Usage
-----
python chembl_validation_spearman.py \
  --rankings path/to/your_drug_rankings.csv \
  --ic50 path/to/SuppData8.xlsx path/to/SuppData3.xlsx \
  --cell-lines MCF7 HepG2 HT29 \
  --score-column score \
  --compound-column compound \
  --join-key compound \
  --out-dir results/

References
----------
- Bin Chen et al., Nature Communications (2017): validation correlates reversal potency with IC50
  in MCF7 (BRCA), HepG2 (LIHC), HT29 (COAD) using Spearman correlation and median IC50. 
"""

import os
import sys
HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))  # go up 3 levels
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from conf import DISEASE, CS_OUT, DATA_DIR, CS_DIR,IMG_DIR,\
    cell_line, diseases_of, LOGS_DIR,\
        cs_filename, disease_run_name, cell_line_run_name,\
    cs_on_LM, cs_mith, selected_cs_run_id, \
    cs_log_filename, lincs_metadata_path, chembl_val_log_filename, ic50_file,\
    IC50_ONLY,CS_TH, IC50_EFF_TH, DRUG_COLLAPSE_METHOD

from logger import append_run_metadata

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr, linregress, ttest_ind

def collapse_cs_profiles_to_drug(cs_df,score_col="connectivity_score",drug_col="pert_id",how="best",):
    """
    Collapse multiple LINCS perturbagens per drug into one row per drug 
    Parameters
    ----------
    cs_df : pd.DataFrame
        Input dataframe.
    score_col : str
        Column to aggregate.
    drug_col : str
        Column identifying the drug.
    how : {'median', 'mean', 'best'}
        Aggregation method:
        - 'median': median score per drug
        - 'mean': mean score per drug
        - 'best': minimum score per drug

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per drug and the same score column name.
    """
    df = cs_df.copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[drug_col, score_col]).copy()

    if how == "median":
        out = df.groupby(drug_col, as_index=False).agg({
            score_col: "median",
        })
    elif how == "mean":
        out = df.groupby(drug_col, as_index=False).agg({
            score_col: "mean",
        })
    elif how == "best":
        out = df.groupby(drug_col, as_index=False).agg({
            score_col: "min",
        })
    elif how == None:
        out=df
    else:
        raise ValueError("how must be one of: 'median', 'mean', 'best'")

    return out

def resolve_cs_run_id(  cs_runs_tsv,  disease_run_id,   drug_run_id,   cs_on_LM,  mith,  selected_cs_run_id=None,):
    """
    Retrieve the cs_run_id from cs log file,
    or return selected_cs_run_id, if provided
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

    hit = runs[
        (runs["disease_run_id"] == disease_run_id) &
        (runs["drug_run_id"] == drug_run_id) &
        (runs["cs_on_LM"] == int(cs_on_LM)) &
        (runs["mith"] == int(mith))
    ].copy()

    if len(hit) == 0:
        raise ValueError(
            "No matching CS run found in cs_runs.tsv for:\n"
            f"disease_run_id={disease_run_id}, "
            f"drug_run_id={drug_run_id}, "
            f"cs_on_LM={int(cs_on_LM)}, "
            f"mith={int(mith)}"
        )

    if len(hit) > 1:
        hit["timestamp"] = pd.to_datetime(hit["timestamp"])
        hit = hit.sort_values("timestamp", ascending=False)

    return hit.iloc[0]["cs_run_id"]

def add_pert_id_to_cs( lincs_metadata_path,cs_df, cs_id_col='LINCS_id', metadata_id_col='id', metadata_pert_col='pert_id'):
    """
    Adds a 'pert_id' column to a CS dataframe by mapping LINCS_id via LINCS metadata.
    
    Parameters
    ----------
    cs_df : pd.DataFrame
        Connectivity score dataframe containing LINCS_id column
    lincs_metadata_path : str
        Path to lincs_sig_info_new.csv
    cs_id_col : str
        Column in cs_df (default 'LINCS_id')
    metadata_id_col : str
        Column in metadata corresponding to LINCS_id (default 'id')
    metadata_pert_col : str
        Column in metadata for pert_id (default 'pert_id')
    """

    if not metadata_pert_col in cs_df.columns:
        # load only what we need
        meta = pd.read_csv(lincs_metadata_path, usecols=[metadata_id_col, metadata_pert_col], dtype='str')
    
        # build mapping dict
        id_to_pert = dict(zip(meta[metadata_id_col], meta[metadata_pert_col]))
    
        # map
        cs_df[metadata_pert_col] = cs_df[cs_id_col].map(id_to_pert)
    
        # optional sanity check
        n_missing = cs_df['pert_id'].isna().sum()
        if n_missing > 0:
            print(f"Warning: {n_missing} LINCS_id values could not be mapped to pert_id")

    return cs_df

    
def load_drug_rankings(path, filename=None, pert_time='all',mith='mith', lincs_metadata_path =None):
    
    
    if not filename:
        filename=mith+'_connectivity_score.tsv'
    df= pd.read_csv(path/filename, sep='\t', header=0, dtype='str')
    
    if not pert_time == 'all':
        return dr[dr.perturbation_time==pert_time]
    
    if lincs_metadata_path:
        df = add_pert_id_to_cs(lincs_metadata_path, df) 
        
        
    return df

def translate_cl(cell_line):
    if cell_line == 'HT29':
        return 'HT-29'
    return cell_line

def load_IC50(ic50_file, cancer_type, cell_line, IC50_ONLY=True):#, median_IC50=False):
    
    cell_line = translate_cl(cell_line)
    df=pd.read_excel(ic50_file,\
                     sheet_name=cancer_type, header =0, usecols=['pert_iname', 'pert_id','standard_value', 'standard_units',\
                                 'standard_type', 'cell_line','activity', 'standard_value_median'])
    #filter for cell line
    df = df[df.cell_line.str.upper()==cell_line]
    
    if IC50_ONLY:
        df=df[df.standard_type=='IC50']
    
    # rename for later merge
    df.rename( columns = {'pert_iname':'drug'}, inplace=True)
    return df.sort_values(by='standard_value_median')

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
    ic50_col="standard_value_median",
    log_ic50_col="log10_ic50",
    drug_label_col="drug",
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
    ax.text(0.03,0.97, txt, transform=ax.transAxes, va="top",ha="left",  bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),    )
    
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

#%%
if __name__=="__main__":
    
    cs_drug_colname = 'pert_id'  #"drug"
    ic50_drug_colname = "pert_id" # "drug"
    
    print('running correlations between chembl IC50 and drug rankings for disease:', DISEASE)
    print(DISEASE, CS_OUT,  cell_line)
    
    # get CS filename for current run
    cs_drug_file= cs_filename
    cs_run_id = resolve_cs_run_id(cs_runs_tsv=cs_log_filename, disease_run_id=disease_run_name,\
    drug_run_id=cell_line_run_name, cs_on_LM=cs_on_LM,   mith=cs_mith,\
        selected_cs_run_id=selected_cs_run_id)

        
    print("Resolved cs_run_id:", cs_run_id)
    cs_drug_file = f"{cs_run_id}.tsv"
    print("Using CS file:", CS_OUT / cs_drug_file)
    
    dr_uncollapsed=load_drug_rankings(CS_OUT, filename = cs_drug_file)
    print(dr_uncollapsed.shape, 'drugs ')
    DRUG_COLLAPSE_METHOD=None
    dr=collapse_cs_profiles_to_drug(cs_df=dr_uncollapsed, drug_col=cs_drug_colname, how=DRUG_COLLAPSE_METHOD)
    print(dr.shape, 'drugs ')

    #%%
    
        #filter for 1 nM
    #%%
    
    # median_IC50=True
    
    ic50=load_IC50(ic50_file, DISEASE, cell_line, IC50_ONLY=IC50_ONLY)#, median_IC50=median_IC50)
    print(ic50.shape, 'drugs with IC50 value for cell line', cell_line)
    
    # calc median ic50
    
    
    #%% merge on common drugs
    merged = pd.merge(dr, ic50, left_on=cs_drug_colname, right_on=ic50_drug_colname, how="inner", suffixes=("_dr","_ic50"))
    
    merged = merged.dropna(subset=["connectivity_score","standard_value_median"]).copy()
    print('merged data on common drugs:', merged.shape)
    
    #%% Optional: check overlap between our overlap of ic50 vs CS score
    # and BinChen's overlap of ic50vs sRGES:
    
    BC_merged = df=pd.read_excel(DATA_DIR/'BinChen2017'/'SD5.xlsx',\
                                 sheet_name=DISEASE)
    overlap_BC = set(merged.drug).intersection(set(BC_merged.pert_iname))
    print(overlap_BC, len(overlap_BC))
    
    #%%# calculate spearman coefficient between
    #  cell_line cells drug rankings vs IC50s
    # Following Bin Chen 2017, more negative reversal score 
    # (RGES/CS) should be associated with stronger efficacy (lower IC50).
    # We report rho for (CS vs IC50) and (CS vs log10(IC50)).
    merged["log10_ic50"] = np.log10(merged["standard_value_median"])
    
    rho_linear, p_linear = spearmanr(merged["connectivity_score"], merged["standard_value_median"])
    rho_log, p_log = spearmanr(merged["connectivity_score"], merged["log10_ic50"])
    
    print('linear IC50: rho=', np.round(rho_linear,2),' pval =' ,np.round(p_linear, 2))
    print('log IC50: rho=', np.round(rho_log, 2),' pval =', np.round(p_log,2))
    
    #%% Optional: check rho with SD5 data:
    BC_merged["log10_ic50"] = np.log10(BC_merged["standard_value"])
    
    rho_linear_bc, p_linear_bc = spearmanr(BC_merged["sRGES"], BC_merged["standard_value"])
    rho_log_bc, p_log_bc = spearmanr(BC_merged["sRGES"], BC_merged["log10_ic50"])
    
    print('linear IC50 SD5: rho=', np.round(rho_linear_bc,2),' pval =' ,np.round(p_linear_bc, 2))
    print('log IC50 SD5: rho=', np.round(rho_log_bc, 2),' pval =', np.round(p_log_bc,2))
    
    #%% Prec/ Rec
    
    classified_df, pr_metrics = classify_ic50_vs_cs(
        merged,
        score_col="connectivity_score",
        ic50_col="standard_value_median",
        cs_threshold=CS_TH,
        ic50_threshold=IC50_EFF_TH,
    )
    print(
        "classification:",
        f"TP={pr_metrics['tp']}, FP={pr_metrics['fp']}, "
        f"FN={pr_metrics['fn']}, TN={pr_metrics['tn']}, "
        f"precision={np.round(pr_metrics['precision'], 2)}, "
        f"recall={np.round(pr_metrics['recall'], 2)}",
    )
    #%% plot

    plot_file = IMG_DIR / f"{DISEASE}_{cs_run_id}_fig3_style.png"
    fig, axes, stats_dict = plot_binchen_fig3_style(
        classified_df = classified_df,
        metrics=pr_metrics,
        disease=DISEASE,
        cell_line_name=cell_line,
        score_col="connectivity_score",
        drug_label_col="drug",
        annotate_top_n=5,
        output_file=plot_file,
        cs_threshold=CS_TH,
        ic50_threshold=IC50_EFF_TH,
    )
#%% Log metadata

    
    run_metadata_data = {
        # identity
        "correlation_run_id": datetime.now().strftime("%d_%m_%Y_%H_%M_%S"),
        "datetime": datetime.now().isoformat(),

        # inputs
        "cs_file": str(CS_OUT / cs_drug_file),
        "ic50_file": str(ic50_file),
        "lincs_metadata_file": str(lincs_metadata_path),
        "disease_run_id": disease_run_name,
        "drug_run_id" : cell_line_run_name,

        # parameters
        "cs_drug_colname": cs_drug_colname,
        "ic50_drug_colname": ic50_drug_colname,
        # "cs_value_col": cs_value_col,
        # "ic50_value_col": ic50_value_col,
        
        # sizes
        "drug_summarization_method": DRUG_COLLAPSE_METHOD,
        "n_cs_rows_before_drug_collapse": len(dr_uncollapsed),
        "n_cs_rows_after_drug_collapse": len(dr),
        "n_ic50_rows": len(ic50),
        "n_merged_rows": len(merged),
        "n_unique_drugs": merged[cs_drug_colname].nunique(),
        
        

        # results
        "spearman_r": np.round(rho_linear,2),
        "spearman_pval": np.round(p_linear,2),
        "spearman_r_log": np.round(rho_log,2),
        "spearman_pval_log": np.round(p_log,2),
        
        # optional: correlations with BC data
        "spearman_r_SD5": np.round(rho_linear_bc,2),
        "spearman_pval_SD5": np.round(p_linear_bc,2),
        "spearman_r_log_SD5": np.round(rho_log_bc,2),
        "spearman_pval_log_SD5": np.round(p_log_bc,2),
        
        # precision recall metrics
        "tp": pr_metrics["tp"],
        "fp": pr_metrics["fp"],
        "fn": pr_metrics["fn"],
        "tn": pr_metrics["tn"],
        "precision": np.round(pr_metrics["precision"], 4),
        "recall": np.round(pr_metrics["recall"], 4),

        # output
        "output_plot_file": str(plot_file)
    }
    

    append_run_metadata(chembl_val_log_filename, run_metadata_data)
    

#%%   
    # #%%
    # # todo: raccogliere i calcoli di sopra in un dizionario che
    # # abbia tutte le voci, e riempirlo per tutte le disease e eprt time
    # # renderlo poi un dizionario, e traslare tutto su ntobeook (lo stesso di cehmbl validation)
    # # per printare il dataframe per bene. 
    
    # # salvare eventualmente in un dizionario a parte,
    # # i valori di IC50_drugs, common_drugs, overlap Drugs
    # data_of = {}
    # for cell_line in ['MCF7','HEPG2','HT29']:
    #     print('-----------------------',cell_line)
    #     disease = diseases_of[cell_line]
    #     ic50 = load_IC50(disease, cell_line)
    #     BC_merged = df=pd.read_excel(DATA_DIR/'BinChen2017'/'SD5.xlsx',\
    #                                  sheet_name=disease)
    #     print(ic50.shape, 'drugs with IC50 value for cell line', cell_line)
    
    #     for pert_time in ['6h', '24h']:
    #         cs_out=CS_DIR/'output'/disease+'_2025_'+pert_time
    #         print(cell_line, 'pert time',pert_time)
    #         dr=load_drug_rankings(cs_out)
    
    #         merged = pd.merge(dr, ic50, on="drug", how="inner", suffixes=("_dr","_ic50"))
    #         merged = merged.dropna(subset=["connectivity_score","standard_value_median"]).copy()
    #         print('merged data on common drugs:', merged.shape)
            
    #         overlap_BC = set(merged.drug).intersection(set(BC_merged.pert_iname))
    #         print('overlap between merged data with MCS vs merged data with sRGES from BinCHen2017', len(overlap_BC))
            
    #         merged["log10_ic50"] = np.log10(merged["standard_value_median"])
    
    #         rho_linear, p_linear = spearmanr(merged["connectivity_score"], merged["standard_value_median"])
    #         rho_log, p_log = spearmanr(merged["connectivity_score"], merged["log10_ic50"])
    
    #         print('linear IC50: rho=', np.round(rho_linear,2),' pval =' ,np.round(p_linear, 2))
    #         print('log IC50: rho=', np.round(rho_log, 2),' pval =', np.round(p_log,2))
            
    #         BC_merged["log10_ic50"] = np.log10(BC_merged["standard_value"])
            
    #         BC_rho_linear, BC_p_linear = spearmanr(BC_merged["sRGES"], BC_merged["standard_value"])
    #         BC_rho_log, BC_p_log = spearmanr(BC_merged["sRGES"], BC_merged["log10_ic50"])
            
    #         print('linear IC50: rho=', np.round(rho_linear,2),' pval =' ,np.round(p_linear, 2))
    #         print('log IC50: rho=', np.round(rho_log, 2),' pval =', np.round(p_log,2))
            
    #         data_of[disease+'_'+pert_time] = [ic50.shape[0], dr.shape[0],\
    #                                           merged.shape[0],BC_merged.shape[0],\
    #                                         len(overlap_BC),rho_linear, p_linear,\
    #                                             rho_log, p_log, BC_rho_linear,\
    #                                             BC_p_linear, BC_rho_log, BC_p_log]
                
    
    # data = pd.DataFrame(data_of).transpose()
    # data.columns= ['n_IC50','n_MCS', 'n_MCSvIC50', \
    #                                         'n_sRGESvIC50','n_overlapSRGESvMCS',\
    #                                         'rho','p','rho_log','p_log',\
    #                                         'rho_sRGES', 'p_sRGES','rho_log_sRGES',\
    #                                             'p_log_sRGES']
