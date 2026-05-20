#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 20:37:14 2025


@ author: L-F-S
calculate RGES connectivity score between disease signature and drug sinatures.


"""
import os
import sys
import numpy as np
import pandas as pd
import time
from scipy import stats
from conf import DISEASE
from loader import load_disease_signature, load_single_drug_signature


# FUNCTIONS


def get_common_genes(disease_signature, drug_signature, id_col='gene_id'):
    '''
    Identifies common items between drug and disease signature.
    Returns indexes of common items in both signatures.
    '''


    missing_items = set(disease_signature[id_col]) - set(drug_signature[id_col])
    addional_drug_items = set(drug_signature[id_col]) - set(disease_signature[id_col])
    common_items = set(disease_signature[id_col]) & set(drug_signature[id_col])


    if not len(common_items) + len(addional_drug_items) == drug_signature.shape[0]:
        raise ValueError(f'common items + additional drug items != len(drug_signature) for {id_col}')
    if not len(common_items) + len(missing_items) == disease_signature.shape[0]:
        raise ValueError(f'common items + missing drug items != len(disease_signature) for {id_col}')


    drug_common = drug_signature[drug_signature[id_col].isin(common_items)].sort_values(by=id_col)
    disease_common = disease_signature[disease_signature[id_col].isin(common_items)].sort_values(by=id_col)


    left_ids = disease_signature.loc[disease_common.index].reset_index()[id_col]
    right_ids = drug_signature.loc[drug_common.index].reset_index()[id_col]
    if not np.all(left_ids.to_numpy() == right_ids.to_numpy()):
        raise ValueError(f'common item indexes actually different between drug and disease for {id_col}!')


    return disease_common.index, drug_common.index


def get_common_genes_old(disease_signature, drug_signature):
    '''
    Identifies common genes between drug and disease signature. Returns indexes
    of common genes in both signatures.
    Input:
        - disease_signature: pd.DataFrame with columns ['gene_id', ...]
        - drug_signature: pd.DataFrame with columns ['gene_id', ...]
    Output:
        - Tuple of DataFrame.index of common genes in disease and drug signatures

    '''


    missing_genes=set.difference(set(disease_signature.gene_id),set(drug_signature.gene_id))
    # print('missing genes from mithril', len(missing_genes))
    addional_drug_items=set.difference(set(drug_signature.gene_id),set(disease_signature.gene_id))
    # print('addional_mith_items', len(addional_drug_items))


    common_genes=set.intersection(set(disease_signature.gene_id),set(drug_signature.gene_id))
    # print('genes in common', len(common_genes))
    if not len(common_genes)+len(addional_drug_items)==drug_signature.shape[0]:
        raise ValueError('common genes+ additional mithril items != len(mithril_gene)')
    if not len(common_genes)+len(missing_genes)==disease_signature.shape[0]:
        raise ValueError('riprova: common genes+ missing mithril items != len(DEG_genes)')


    # Get subset of signatures with common genes
    drug_common_genes_signatures=drug_signature[drug_signature['gene_id'].isin(common_genes)].sort_values(by='gene_id')
    disease_common_genes_signatures=disease_signature[disease_signature.gene_id.isin(common_genes)].sort_values(by='gene_id')

    if not np.unique(disease_signature.loc[disease_common_genes_signatures.index].reset_index()['gene_id']==drug_signature.loc[drug_common_genes_signatures.index].reset_index()['gene_id']):
        raise ValueError('common gene indexes actually different between drug and disease!')
    return disease_common_genes_signatures.index, drug_common_genes_signatures.index



def rank_genes(drug_signature, disease_signature, drug_col_name, disease_col_name, id_col='gene_id'):
    '''
    Compute the array V of sorted drug item indexes, indexed by sorted disease
    indexes.
    Input:
        - drug_signature: pd.DataFrame with columns [id_col, drug_col_name, ...]
        - disease_signature: pd.DataFrame with columns [id_col, disease_col_name, ...]
        - drug_col_name: str, column name to sort drug_signature by
        - disease_col_name: str, column name to sort disease_signature by
        - id_col: str, column name of gene IDs
     Output:
        NumPy array where the indices are sorted disease gene indexes
        and the values are indexes of corresponding sorted drug genes.
    '''
    sorted_drug_genes = drug_signature.sort_values(by=drug_col_name, ascending=False).reset_index(drop=True)
    sorted_disease_genes = disease_signature.sort_values(by=disease_col_name)


    merged_data_on_disease_indexes = pd.merge(sorted_disease_genes,sorted_drug_genes.reset_index(), on=id_col)


    V = merged_data_on_disease_indexes['index'].to_numpy()
    return V


def compute_KS_disease_sorted(disease_disregulated_genes, V, drug_genes):
    '''
    Compute Kolmogorov-Smirnov (KS) like statistic.
    From Lamb et al., 2006 supplementary
    Input:
       - disease_disregulated_genes: list or array of gene IDs
       - V: NumPy array of drug gene indices sorted on disease indices
       - drug_genes: list or array of gene IDs in the drug signature
    Output:
       - a: float, maximum positive difference
       - b: float, maximum negative difference
       - s: int, number of disregulated genes
       - r: int, total number of genes in reference drug expression data
    '''
    # number of disregulated genes:
    s=len(disease_disregulated_genes)


    # total number of genes in reference drug expression data:
    r=len(drug_genes)

    V_over_r = V/r
    a = np.max(np.arange(1, s+1)/s - V_over_r)
    b = np.max(V_over_r - (np.arange(s)/s))

    return a, b, s, r


def compute_KS(disease_disregulated_genes, V, drug_genes):
    '''
    Compute Kolmogorov-Smirnov (KS) statistic.
    From Lamb et al., 2006 supplementary
    Input:
       - disease_disregulated_genes: list or array of gene IDs
       - V: NumPy array of sorted drug gene indices
       - drug_genes: list or array of gene IDs in the drug signature
    Output:
       - a: float, maximum positive difference
       - b: float, maximum negative difference
       - s: int, number of disregulated genes
       - r: int, total number of genes in reference drug expression data
    '''
    V=np.sort(V)
    # number of disregulated genes:
    s=len(disease_disregulated_genes)


    # total number of genes in reference drug expression data:
    r=len(drug_genes)

    V_over_r = V/r
    a = np.max(np.arange(1, s+1)/s - V_over_r)
    b = np.max(V_over_r - (np.arange(s)/s))

    return a, b, s, r




def compute_KS_evil_twin(disease_disregulated_genes, V, drug_genes):
    '''
    Compute Kolmogorov-Smirnov (KS) statistic.
    From Lamb et al., 2006 supplementary
    Input:
       - disease_disregulated_genes: list or array of gene IDs
       - V: NumPy array of sorted drug gene indices
       - drug_genes: list or array of gene IDs in the drug signature
    Output:
       - a: float, maximum positive difference
       - b: float, maximum negative difference
       - s: int, number of disregulated genes
       - r: int, total number of genes in reference drug expression data
    '''

    # number of disregulated genes:
    s=len(disease_disregulated_genes)


    # total number of genes in reference drug expression data:
    r=len(drug_genes)

    V_over_r = V/r
    a = np.max(np.abs(np.arange(1, s+1)/s - V_over_r))
    b = np.max(np.abs(V_over_r - (np.arange(s)/s)))
    ks=np.max(a,b)

    return ks, s, r




def lamb_normalize(cs_list):
    '''
    CS normalization following Lamb et al., 2006
    Input:
        cs_list: list, connectivity scores for drugs of interest
    Output:
        list, normalized connectivity scores
    '''

    p=max(cs_list)
    q=min(cs_list)

    def extr(x):
        if x>0:
            return p
        return -q

    return [x/extr(x) for x in cs_list]




def calculate_CS(a_up, a_down, b_up, b_down):
    '''
    Calculate connectivity score.
    from Lamb et al., 2006
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - CS: float, CS value in the interval [-1, 1]
    '''

    if a_up>b_up:
        ks_up=a_up
    else:
        ks_up=-b_up

    if a_down>b_down:
        ks_down=a_down
    else:
        ks_down=-b_down

    if ks_up*ks_down>0:
        return 0
    return ks_up-ks_down




def calculate_CS_evil_twin(ks_up, ks_down):
    '''
    Calculate connectivity score. But its evil twin.
    For testing
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - CS: float, CS value in the interval [-1, 1]
    '''
    return ks_up-ks_down




def calculate_RGES(a_up, a_down, b_up, b_down):
    '''
    Calculate Reverse Gene Expression Score (RGES).
    from Bin Chen et al., 2017
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - RGES: float, RGES value in the interval [-2, 2]
    '''

    if a_up>b_up:
        ks_up=a_up
    else:
        ks_up=-b_up

    if a_down>b_down:
        ks_down=a_down
    else:
        ks_down=-b_down

    return ks_up-ks_down


def calculate_pval_weighted_sRGES(rges_values, p_values, p_cutoff=0.05):
    '''
    NetCoS-specific p-value-weighted summarised RGES (NOT the original Bin Chen
    2017 sRGES).  Kept for backward compatibility with earlier NetCoS
    experiments.  See `calculate_bin_chen_sRGES` for the published formula.

    Definition (NetCoS variant):
        sRGES = weighted mean of per-perturbation RGES values, restricted to
        instances that significantly reverse the disease signature:
            - RGES < 0  (drug reverses disease direction)
            - p_value < p_cutoff
        Each qualifying instance is weighted by -log10(p_value).

    Input:
        - rges_values: array-like of float, per-perturbation-instance RGES values
        - p_values: array-like of float, per-perturbation-instance p-values from
                    Monte Carlo simulation
        - p_cutoff: float, significance threshold (default: 0.05)
    Output:
        - sRGES: float, summarised RGES in (-inf, 0]; 0.0 if no
                 significant reversals exist across all instances of the drug
    '''
    rges  = np.asarray(rges_values,  dtype=float)
    pvals = np.asarray(p_values, dtype=float)

    mask = (rges < 0) & (pvals < p_cutoff)
    if not np.any(mask):
        return 0.0

    weights = -np.log10(pvals[mask])
    return float(np.average(rges[mask], weights=weights))


# Backward-compatibility alias for the older (mis-named) function.
# Earlier code referred to the p-value-weighted variant above as
# `calculate_sRGES`.  The published Bin Chen 2017 sRGES is implemented as
# `calculate_bin_chen_sRGES` further below.
calculate_sRGES = calculate_pval_weighted_sRGES


def _bin_chen_diff(df, drug_id_col, cell_col, dose_col, time_col, cs_col,
                   ref_dose=10.0, ref_time=24, time_values=(6, 24),
                   dose_min=0.0):
    '''
    Compute the dose-time correction vector `diff` used by Bin Chen's
    sRGES (Bin Chen et al., 2017; sRGES_all_cmpds.R / core_functions.R).

    The procedure replicates the R reference implementation:
        1) restrict to profiles with pert_dose > dose_min and
           pert_time in time_values (default 6h and 24h)
        2) for every (drug, cell) self-pair where the reference profile
           is at ref_dose (10 microM) and ref_time (24h), compute
                cmap_diff = RGES_reference - RGES_test
        3) average cmap_diff per (dose_bin, time) cell, where
           dose_bin = "low"  if pert_dose < ref_dose
                    = "high" otherwise.

    Input:
        df              : pd.DataFrame, one row per LINCS perturbation instance
        drug_id_col     : column name of the drug identifier (e.g. pert_iname)
        cell_col        : column name of the cell-line identifier
        dose_col        : column name of the perturbation dose
        time_col        : column name of the perturbation time
        cs_col          : column name of the per-instance RGES / connectivity
                          score (cmap_score in the BinChen LINCS file)
        ref_dose        : reference dose, default 10.0 microM
        ref_time        : reference time, default 24 h
        time_values     : tuple of allowed times after filtering, default (6, 24)
        dose_min        : strict lower bound on pert_dose, default 0.0

    Output:
        dict with keys ("high",24), ("high",6), ("low",24), ("low",6) ->
        mean cmap_diff (0.0 for missing cells).  Cells not represented in the
        data default to 0.0 so the adjustment becomes a no-op.
    '''
    subset = df[(df[dose_col] > dose_min) & df[time_col].isin(time_values)].copy()
    if subset.empty:
        return {('high', ref_time): 0.0, ('high', 6): 0.0,
                ('low',  ref_time): 0.0, ('low',  6): 0.0}

    # use a row-identifier so we can drop self-self pairs (id.x != id.y)
    subset = subset.reset_index(drop=False).rename(columns={'index': '_row'})
    ref = subset[(subset[time_col] == ref_time) & (subset[dose_col] == ref_dose)]

    if ref.empty:
        return {('high', ref_time): 0.0, ('high', 6): 0.0,
                ('low',  ref_time): 0.0, ('low',  6): 0.0}

    # self-merge on (drug, cell) — reference (.x) vs all other profiles (.y)
    pairs = ref.merge(subset, on=[drug_id_col, cell_col], suffixes=('_x', '_y'))
    pairs = pairs[pairs['_row_x'] != pairs['_row_y']]
    if pairs.empty:
        return {('high', ref_time): 0.0, ('high', 6): 0.0,
                ('low',  ref_time): 0.0, ('low',  6): 0.0}

    pairs['cmap_diff'] = pairs[f'{cs_col}_x'] - pairs[f'{cs_col}_y']
    pairs['dose_bin']  = np.where(pairs[f'{dose_col}_y'] < ref_dose, 'low', 'high')
    pairs['time_y']    = pairs[f'{time_col}_y'].astype(int)

    diff = pairs.groupby(['dose_bin', 'time_y'])['cmap_diff'].mean().to_dict()

    # ensure every key exists (default 0.0)
    for key in [('high', ref_time), ('high', 6), ('low', ref_time), ('low', 6)]:
        diff.setdefault(key, 0.0)
    return diff


def _apply_bin_chen_getsRGES(rges, cor, pert_dose, pert_time, diff, max_cor,
                             ref_dose=10.0, ref_time=24,
                             match_r_indexing=True):
    '''
    Per-instance Bin Chen 2017 RGES adjustment (Python port of getsRGES in
    Bin-Chen-Lab/RGES/core_functions.R).

    Input:
        rges            : float, raw per-instance RGES (cmap_score)
        cor             : float, cell line correlation with target tumour (TCGA)
        pert_dose       : float, perturbation dose
        pert_time       : float, perturbation time
        diff            : dict, output of `_bin_chen_diff`
        max_cor         : float, max cor across all instances of the disease
        ref_dose        : float, dose threshold for "low"/"high" binning (10)
        ref_time        : float, time threshold for "short"/"long" binning (24)
        match_r_indexing: bool, default True. Replicate the reference R code
                          exactly, including the apparent index mismatch in
                          core_functions.R `getsRGES` (see comment below). Set
                          to False to use a "fixed" semantically-correct
                          indexing.

    The reference R `getsRGES` (line-for-line):
        if (pert_time == "short" & pert_dose == "low")  -> sRGES += diff[4]
        if (pert_dose == "low"  & pert_time == "long")  -> sRGES += diff[2]
        if (pert_dose == "high" & pert_time == "short") -> sRGES += diff[1]
        # high/long (reference) -> no adjustment
        return sRGES * cor / max_cor

    In R, `diff` is built via `tapply(..., paste(dose_bin, pert_time), mean)`
    and its names sort alphabetically:
        diff[1] = "high 24"   diff[2] = "high 6"
        diff[3] = "low 24"    diff[4] = "low 6"
    With `match_r_indexing=True` we apply exactly the same combinations as the
    R code (so we reproduce the reference `lincs_cancer_sRGES.csv` bit-for-bit).
    '''
    pert_time_bin = 'short' if pert_time < ref_time else 'long'
    pert_dose_bin = 'low'   if pert_dose < ref_dose else 'high'

    srges = float(rges)
    if match_r_indexing:
        # exact port of the published R `getsRGES`
        if pert_time_bin == 'short' and pert_dose_bin == 'low':
            srges = srges + diff.get(('low',  6),        0.0)  # R diff[4]
        if pert_dose_bin == 'low'   and pert_time_bin == 'long':
            srges = srges + diff.get(('high', 6),        0.0)  # R diff[2]
        if pert_dose_bin == 'high'  and pert_time_bin == 'short':
            srges = srges + diff.get(('high', ref_time), 0.0)  # R diff[1]
    else:
        # semantically-aligned variant (each bin gets its own diff)
        if pert_time_bin == 'short' and pert_dose_bin == 'low':
            srges = srges + diff.get(('low',  6),        0.0)
        if pert_dose_bin == 'low'   and pert_time_bin == 'long':
            srges = srges + diff.get(('low',  ref_time), 0.0)
        if pert_dose_bin == 'high'  and pert_time_bin == 'short':
            srges = srges + diff.get(('high', 6),        0.0)
        # high/long is the reference - no adjustment

    if max_cor and not np.isnan(max_cor):
        return srges * cor / max_cor
    return srges


def calculate_bin_chen_sRGES(df,
                             drug_id_col='pert_iname',
                             cs_col='cmap_score',
                             cell_col='cell_id',
                             dose_col='pert_dose',
                             time_col='pert_time',
                             cor_col='cor',
                             cell_cor_map=None,
                             ref_dose=10.0,
                             ref_time=24,
                             time_values=(6, 24),
                             dose_min=0.0,
                             match_r_indexing=True,
                             return_adjusted=False):
    '''
    Per-drug Summarised RGES (sRGES) following Bin Chen et al., 2017
    (https://github.com/Bin-Chen-Lab/RGES, sRGES_all_cmpds.R + core_functions.R).

    Pipeline:
        1) restrict to profiles with pert_dose > dose_min and
           pert_time in time_values (default {6, 24}h)
        2) merge the per-instance cell-line correlation (cor) — either from
           an existing column `cor_col`, or from the dict `cell_cor_map`
           (cell_id -> cor); rows with missing cor are dropped (this matches
           the R inner-merge to cell_line_<cancer>_tacle.csv)
        3) compute the dose-time correction vector `diff` from same-(drug,
           cell) reference-vs-test pairs
        4) for every instance, adjust RGES via `getsRGES`:
                 sRGES_i = (RGES_i + diff[dose_bin, time_bin]) * cor_i / max(cor)
        5) per drug, sRGES = mean(sRGES_i)

    Input:
        df              : pd.DataFrame, one row per LINCS perturbation instance.
                          Must contain at least drug_id_col, cs_col, cell_col,
                          dose_col, time_col.  cor is optional (see below).
        drug_id_col     : drug identifier column (default 'pert_iname')
        cs_col          : per-instance RGES column (default 'cmap_score')
        cell_col        : cell-line column (default 'cell_id')
        dose_col        : dose column (default 'pert_dose')
        time_col        : time column (default 'pert_time')
        cor_col         : column with TCGA-cell-line correlation (default 'cor').
                          If absent, falls back to `cell_cor_map`, else to
                          cor=1 for every row (no cor weighting).
        cell_cor_map    : optional dict {cell_id -> cor}; used if cor_col not
                          in df.
        ref_dose        : reference dose, default 10.0 microM
        ref_time        : reference time, default 24 h
        time_values     : tuple of admissible times, default (6, 24)
        dose_min        : strict lower bound on pert_dose, default 0.0
        match_r_indexing: bool, default True (reproduce reference R code
                          including its diff-index quirk; required to match
                          BinChen2017 reference output).
        return_adjusted : bool, default False.  If True also return the
                          per-instance adjusted-RGES DataFrame (useful for
                          debugging or downstream weighting).

    Output:
        pd.DataFrame indexed implicitly with columns:
            [drug_id_col, 'mean', 'n', 'median', 'sd', 'sRGES']
        — column order matches BinChen2017 lincs_cancer_sRGES.csv,
        with sRGES == mean of the per-instance adjusted RGES.

        If return_adjusted=True, returns a tuple (per_drug_df, per_instance_df)
        where per_instance_df has the original columns plus 'adjusted_RGES'.
    '''
    if df.empty:
        cols = [drug_id_col, 'mean', 'n', 'median', 'sd', 'sRGES']
        empty = pd.DataFrame(columns=cols)
        if return_adjusted:
            return empty, df.assign(adjusted_RGES=[])
        return empty

    work = df.copy()

    # 1) attach cor
    if cor_col not in work.columns:
        if cell_cor_map is not None:
            work[cor_col] = work[cell_col].map(cell_cor_map)
        else:
            work[cor_col] = 1.0

    # 2) inner-merge semantics: drop rows with missing cor
    work = work.dropna(subset=[cor_col]).copy()
    if work.empty:
        cols = [drug_id_col, 'mean', 'n', 'median', 'sd', 'sRGES']
        empty = pd.DataFrame(columns=cols)
        if return_adjusted:
            return empty, work.assign(adjusted_RGES=[])
        return empty

    # 3) compute diff vector from the FULL (pre-cor-merge) df, mirroring the R
    # code which computes `diff` from `lincs_drug_prediction` directly (before
    # the cor merge).  Use df, not work, for the diff computation.
    diff = _bin_chen_diff(df, drug_id_col=drug_id_col, cell_col=cell_col,
                          dose_col=dose_col, time_col=time_col, cs_col=cs_col,
                          ref_dose=ref_dose, ref_time=ref_time,
                          time_values=time_values, dose_min=dose_min)

    # 4) per-instance adjustment.  max_cor is computed AFTER cor-merge to
    # match R: `max(pred$cor)` where pred is the merged frame.
    max_cor = work[cor_col].max()

    work['adjusted_RGES'] = [
        _apply_bin_chen_getsRGES(row_rges, row_cor, row_dose, row_time,
                                 diff, max_cor,
                                 ref_dose=ref_dose, ref_time=ref_time,
                                 match_r_indexing=match_r_indexing)
        for row_rges, row_cor, row_dose, row_time in zip(
            work[cs_col].to_numpy(dtype=float),
            work[cor_col].to_numpy(dtype=float),
            work[dose_col].to_numpy(dtype=float),
            work[time_col].to_numpy(dtype=float),
        )
    ]

    # 5) per-drug summary; sRGES = mean of adjusted RGES (matches R ddply)
    summary = (work.groupby(drug_id_col)['adjusted_RGES']
                   .agg(mean='mean', n='size', median='median', sd='std')
                   .reset_index())
    summary['sRGES'] = summary['mean']
    summary = summary.sort_values('sRGES').reset_index(drop=True)

    if return_adjusted:
        return summary, work
    return summary


def drug_collapse(cs_df, method='sRGES', drug_id_col='pert_id',
                  cs_col='connectivity_score', p_col='cs_p_value',
                  p_cutoff=0.05,
                  # --- additional kwargs used only by method='sRGES' ----
                  cell_col='cell_id',
                  dose_col='pert_dose',
                  time_col='pert_time',
                  cor_col='cor',
                  cell_cor_map=None,
                  ref_dose=10.0,
                  ref_time=24,
                  match_r_indexing=True):
    '''
    Collapse per-perturbation-instance connectivity scores to a single
    per-drug score.

    Each row of cs_df is one LINCS perturbation instance; this function
    groups by drug identifier and applies the chosen collapse method.

    Input:
        - cs_df: pd.DataFrame with one row per LINCS perturbation instance.
                 Must contain at minimum: drug_id_col, cs_col.  For
                 method='sRGES' must also contain cell_col, dose_col,
                 time_col (cor is optional — see `calculate_bin_chen_sRGES`).
                 For method='pval_sRGES' must also contain p_col.
        - method: str or callable.  Collapse method.  Options:
                  'sRGES'      — Summarised RGES, Bin Chen et al. 2017
                                 (sRGES_all_cmpds.R / core_functions.R).
                                 Dose- and time-corrected per-instance RGES,
                                 averaged per drug. [default]
                  'pval_sRGES' — NetCoS p-value-weighted variant:
                                 weighted mean of significant reversals
                                 (RGES < 0, p < p_cutoff), weighted by
                                 -log10(p).  Earlier in NetCoS this was the
                                 'sRGES' option.
                  'mean'       — arithmetic mean of all instances
                  'median'     — median of all instances
                  'min'        — minimum (most negative) score across all
                                 instances
                  'best'       — alias for 'min'
                  callable     — called as fn(rges_array, pval_array) -> float,
                                 where rges_array and pval_array are 1-D NumPy
                                 float arrays for all instances of a single
                                 drug
        - drug_id_col: str, column to group by (drug identifier).
                       Default 'pert_id'.
        - cs_col: str, connectivity score column. Default 'connectivity_score'.
        - p_col: str, p-value column. Default 'cs_p_value'.
        - p_cutoff: float, significance threshold used by 'pval_sRGES'.
                    Default 0.05.

        sRGES-specific kwargs (ignored for other methods):
        - cell_col, dose_col, time_col, cor_col:
            column names for cell-line, dose, time, and cell-line correlation.
        - cell_cor_map: optional dict {cell_id -> cor}; used if cor_col not
                        in cs_df.
        - ref_dose, ref_time: reference dose/time (10 microM / 24 h).
        - match_r_indexing: replicate the published R `getsRGES` index quirk
                            (default True — required to reproduce BinChen2017
                            reference sRGES files exactly).

    Output:
        For 'sRGES':
            pd.DataFrame with columns [drug_id_col, 'mean', 'n', 'median',
            'sd', 'sRGES', 'collapsed_score', 'n_instances'], where
            'collapsed_score' == 'sRGES' and 'n_instances' == 'n' (so this
            output is drop-in-compatible with the other methods).
        For all other methods:
            pd.DataFrame with columns [drug_id_col, 'collapsed_score',
            'n_instances'], sorted by 'collapsed_score' ascending
            (most negative = best drug candidate).
    '''
    # ---- sRGES (Bin Chen 2017): needs full DataFrame, not per-drug arrays ----
    if method == 'sRGES':
        summary = calculate_bin_chen_sRGES(
            cs_df,
            drug_id_col=drug_id_col,
            cs_col=cs_col,
            cell_col=cell_col,
            dose_col=dose_col,
            time_col=time_col,
            cor_col=cor_col,
            cell_cor_map=cell_cor_map,
            ref_dose=ref_dose,
            ref_time=ref_time,
            match_r_indexing=match_r_indexing,
        )
        # drop-in compatibility columns for callers that expect the legacy
        # (collapsed_score, n_instances) schema
        summary['collapsed_score'] = summary['sRGES']
        summary['n_instances']     = summary['n']
        return summary

    # ---- per-drug-array methods (legacy interface) ----
    if method == 'pval_sRGES':
        collapse_fn = lambda rges, pvals: calculate_pval_weighted_sRGES(
            rges, pvals, p_cutoff=p_cutoff)
    elif method == 'mean':
        collapse_fn = lambda rges, pvals: float(np.mean(rges))
    elif method == 'median':
        collapse_fn = lambda rges, pvals: float(np.median(rges))
    elif method in ('min', 'best'):
        collapse_fn = lambda rges, pvals: float(np.min(rges))
    elif callable(method):
        collapse_fn = method
    else:
        raise ValueError(
            f"Unknown collapse method: {method!r}. "
            f"Options: 'sRGES', 'pval_sRGES', 'mean', 'median', 'min', "
            f"'best', or a callable."
        )

    # p_col is only required for callables / pval_sRGES; allow its absence
    # for the simple summary methods so callers don't have to fabricate one.
    has_p = p_col in cs_df.columns

    records = []
    for drug_id, group in cs_df.groupby(drug_id_col):
        rges  = group[cs_col].to_numpy(dtype=float)
        pvals = (group[p_col].to_numpy(dtype=float)
                 if has_p else np.full(len(group), np.nan))
        score = collapse_fn(rges, pvals)
        records.append({
            drug_id_col: drug_id,
            'collapsed_score': score,
            'n_instances': len(group),
        })

    return pd.DataFrame(records).sort_values('collapsed_score').reset_index(drop=True)


def montecarlo_connectivity(s_up, s_down, r, n_iterations=1000, score_type='bin_chen'):
    '''
    Randomly sample RGES
    Input:
        - s_up: int, number of upregulated genes
        - s_down: int, number of downregulated genes
        - r: int, total number of genes in reference drug expression data
        - n_iterations: int, number of Monte Carlo iterations (default: 1000)
        - score_type: str, optional, which calculation to use.
            options: ['bin_chen', 'sirota', 'lamb'], default: 'bin_chen'
    Output:
        - list of float, sampled RGES values
    '''

    random_RGES_list=[]
    drug_genes=np.arange(r)

    for i in range(n_iterations):

        random_idx = np.random.choice(r, s_up + s_down, replace=False)


        random_V_up = np.sort(random_idx[:s_up])
        random_V_down = np.sort(random_idx[s_up:])

        # Compute random KS stats:
        random_a_up, random_b_up, _, _ = compute_KS(np.arange(s_up), random_V_up, drug_genes)
        random_a_down , random_b_down, _, _ = compute_KS(np.arange(s_down), random_V_down, drug_genes)


        # calculate RGES evil twin:
        if score_type=='evil_twin':
            ks_up,  _, _ = compute_KS(np.arange(s_up), random_V_up, drug_genes)
            ks_down, _, _ = compute_KS(np.arange(s_down), random_V_down, drug_genes)
            random_RGES_list.append(calculate_CS_evil_twin(ks_up, ks_down, random_b_up, random_b_down))

        # Calculate random RGES:
        if score_type=='bin_chen':
            random_RGES_list.append(calculate_RGES(random_a_up, random_a_down, random_b_up, random_b_down))
        if (score_type=='lamb') or(score_type=='sirota') :
            random_RGES_list.append(calculate_CS(random_a_up, random_a_down, random_b_up, random_b_down))


    if score_type=='lamb':
        random_RGES_list=lamb_normalize(random_RGES_list)


    return random_RGES_list




def montecarlo_connectivity_disease_sorted(s_up, s_down, r, n_iterations=1000, score_type='bin_chen'):
    '''
    Randomly sample disease-sorted RGES
    Input:
        - s_up: int, number of upregulated genes
        - s_down: int, number of downregulated genes
        - r: int, total number of genes in reference drug expression data
        - n_iterations: int, number of Monte Carlo iterations (default: 1000)
        - score_type: str, optional, which calculation to use.
            options: ['bin_chen', 'sirota', 'lamb'], default: 'bin_chen'
    Output:
        - list of float, sampled RGES values
    '''

    random_RGES_list=[]
    drug_genes=np.arange(r)

    for i in range(n_iterations):

        # Generate random indexes for up and downregulated genes
        random_up_and_down_indexes=np.random.choice(r, s_up + s_down, replace=False)

        # Create random index dictionary:
        random_V_up=random_up_and_down_indexes[:s_up]
        random_V_down=random_up_and_down_indexes[-s_down:]

        # Compute random KS stats:
        random_a_up, random_b_up, _, _ = compute_KS_disease_sorted(random_up_and_down_indexes[:s_up], random_V_up, drug_genes)
        random_a_down , random_b_down, _, _ = compute_KS_disease_sorted(random_up_and_down_indexes[-s_down:], random_V_down, drug_genes)


        # calculate RGES evil twin:
        if score_type=='evil_twin':
            ks_up,  _, _ = compute_KS_disease_sorted(random_up_and_down_indexes[:s_up], random_V_up, drug_genes)
            ks_down, _, _ = compute_KS_disease_sorted(random_up_and_down_indexes[-s_down:], random_V_down, drug_genes)
            random_RGES_list.append(calculate_CS_evil_twin(ks_up, ks_down, random_b_up, random_b_down))

        # Calculate random RGES:
        if score_type=='bin_chen':
            random_RGES_list.append(calculate_RGES(random_a_up, random_a_down, random_b_up, random_b_down))
        if (score_type=='lamb') or(score_type=='sirota') :
            random_RGES_list.append(calculate_CS(random_a_up, random_a_down, random_b_up, random_b_down))


    if score_type=='lamb':
        random_RGES_list=lamb_normalize(random_RGES_list)


    return random_RGES_list


def random_disease_bin_chen_connectivity(drug_signature, rank_on='magnitude', id_col='gene_id'):
    '''temp
    calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017, of a drug signature versus
   a random genes  ranking
    Input:
        - drug_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - rank_on: str ranking column. default='magnitude'
                options: {'p_value', 'magnitude'}. Do not change unless you
                know what you're doing'

    Output:
        -   measured_RGES: float, calculated RGES value
        -   p_value: float, p-value from Monte Carlo simulation
    '''


    if rank_on=='p_value':
        ranking_col_name=drug_signature.columns[2]
    if rank_on=='magnitude':
        ranking_col_name=drug_signature.columns[1]

    # generate random up and down regulated disease genes
    r= drug_signature.shape[0]
    l=np.random.randint(1,r)
    s_up = r-l
    s_down = l
    gene_id_list = drug_signature[id_col]
    shuffled_gene_ids = gene_id_list[np.random.choice(gene_id_list.index, len(gene_id_list), replace=False)]
    random_disease_up = shuffled_gene_ids[:s_up]
    random_disease_down = shuffled_gene_ids[-s_down:]

    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, pd.DataFrame({id_col: random_disease_up})  , ranking_col_name, id_col, id_col=id_col)
    V_down = rank_genes(drug_signature, pd.DataFrame({id_col: random_disease_down}), ranking_col_name, id_col, id_col=id_col)

    # Compute random KS stats:
    a_up, b_up, _, _ = compute_KS(random_disease_up, V_up, drug_signature[id_col])
    a_down , b_down, _, _ = compute_KS(random_disease_down, V_down, drug_signature[id_col])

    # Compute RGES:
    measured_RGES_for_random_disease = calculate_RGES(a_up, a_down, b_up, b_down)

    return measured_RGES_for_random_disease




def bin_chen_connectivity(disease_signature, drug_signature, rank_on='magnitude', id_col='gene_id',\
    drug_ranking_col_name=None, disease_ranking_col_name=None):
    '''
    calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017
    Input:
        - disease_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - drug_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - rank_on: str ranking column. default='magnitude'
                options: {'p_value', 'magnitude'}. Do not change unless you
                know what you're doing'

    Output:
        -   measured_RGES: float, calculated RGES value
        -   p_value: float, p-value from Monte Carlo simulation
    '''

    if drug_ranking_col_name is None or disease_ranking_col_name is None:
        if rank_on=='p_value':
            drug_ranking_col_name=drug_signature.columns[2]
            disease_ranking_col_name=disease_signature.columns[2]
        if rank_on=='magnitude':
            drug_ranking_col_name=drug_signature.columns[1]
            disease_ranking_col_name=disease_signature.columns[1]

    # # Get lists of up (down) regulated genes: older with 2vs
    disease_signature_up = disease_signature[disease_signature[disease_ranking_col_name]>0]
    disease_signature_down = disease_signature[disease_signature[disease_ranking_col_name]<0]

    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, disease_signature_up, drug_ranking_col_name, disease_ranking_col_name, id_col=id_col)
    V_down = rank_genes(drug_signature, disease_signature_down, drug_ranking_col_name, disease_ranking_col_name, id_col=id_col)

    # Compute KS statistic for up and down regulated genes:
    a_up, b_up, s_up, r = compute_KS(disease_signature_up[id_col], V_up, drug_signature[id_col])
    a_down , b_down, s_down, r = compute_KS(disease_signature_down[id_col], V_down, drug_signature[id_col])

    # Compute RGES:
    measured_RGES = calculate_RGES(a_up, a_down, b_up, b_down)

    # Calculate two tailed p-value (no assumption on the direction of RGES
    # between disease and drug) for measured RGES, using random sampling:
    n_iterations=1000
    random_RGES_list = montecarlo_connectivity(s_up, s_down, r, n_iterations, score_type='bin_chen')
    p_value=np.sum(np.abs(np.array(random_RGES_list))>np.abs(measured_RGES))/n_iterations
    return measured_RGES, p_value


def bin_chen_connectivity_sorted(disease_signature, drug_signature, rank_on='magnitude', id_col='gene_id',\
    drug_ranking_col_name=None,    disease_ranking_col_name=None):
    '''
    calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017
    Input:
        - disease_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - drug_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - rank_on: str ranking column. default='magnitude'
                options: {'p_value', 'magnitude'}. Do not change unless you
                know what you're doing'

    Output:
        -   measured_RGES: float, calculated RGES value
        -   p_value: float, p-value from Monte Carlo simulation
    '''

    if drug_ranking_col_name is None or disease_ranking_col_name is None:
        if rank_on=='p_value':
            drug_ranking_col_name=drug_signature.columns[2]
            disease_ranking_col_name=disease_signature.columns[2]
        if rank_on=='magnitude':
            drug_ranking_col_name=drug_signature.columns[1]
            disease_ranking_col_name=disease_signature.columns[1]

    # # Get lists of up (down) regulated genes: older with 2vs
    disease_signature_up = disease_signature[disease_signature[disease_ranking_col_name]>0]
    disease_signature_down = disease_signature[disease_signature[disease_ranking_col_name]<0]

    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, disease_signature_up, drug_ranking_col_name, disease_ranking_col_name, id_col=id_col)
    V_down = rank_genes(drug_signature, disease_signature_down, drug_ranking_col_name, disease_ranking_col_name, id_col=id_col)

    # Compute KS statistic for up and down regulated genes:
    a_up, b_up, s_up, r = compute_KS_disease_sorted(disease_signature_up[id_col], V_up, drug_signature[id_col])
    a_down , b_down, s_down, r = compute_KS_disease_sorted(disease_signature_down[id_col], V_down, drug_signature[id_col])

    # Compute RGES:
    measured_RGES = calculate_RGES(a_up, a_down, b_up, b_down)

    # Calculate two tailed p-value (no assumption on the direction of RGES
    # between disease and drug) for measured RGES, using random sampling:
    n_iterations=1000
    random_RGES_list = montecarlo_connectivity(s_up, s_down, r, n_iterations, score_type='bin_chen')
    p_value=np.sum(np.abs(np.array(random_RGES_list))>np.abs(measured_RGES))/n_iterations
    return measured_RGES, p_value


calc_connectivity_score_with = {"bin_chen":bin_chen_connectivity, "bin_chen_disease_sorted":bin_chen_connectivity_sorted, "lamb":"TODO: implement me!"}


calc_drug_collapse_with = {
    'sRGES':  lambda cs_df, **kw: drug_collapse(cs_df, method='sRGES',  **kw),
    'mean':   lambda cs_df, **kw: drug_collapse(cs_df, method='mean',   **kw),
    'median': lambda cs_df, **kw: drug_collapse(cs_df, method='median', **kw),
    'min':    lambda cs_df, **kw: drug_collapse(cs_df, method='min',    **kw),
    'best':   lambda cs_df, **kw: drug_collapse(cs_df, method='min',    **kw),
}
