#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Spearman correlation between a drug-ranking score (e.g., RGES/sRGES/NetCos)
and ChEMBL IC50 values for MCF7, HepG2, and HT29 (Bin Chen 2017 validation).

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

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from conf import DISEASE, CS_OUT, DATA_DIR, TSR_OUT_DRUG, TSR_OUT_DISEASE, CS_IN_DISEASE, CS_IN_DRUG
print('running correlations between chembl IC50 and drug rankings for disease:', DISEASE)

# select 
def load_drug_rankings(path, filename=None, pert_time='all',mith='mith'):
    
    
    if not filename:
        filename=mith+'_connectivity_score.tsv'
    df= pd.read_csv(path+filename, sep='\t', header=0, usecols=[0,1,2,3])
    
    if not pert_time == 'all':
        return dr[dr.perturbation_time==pert_time]
    
    return df
    
dr=load_drug_rankings(CS_OUT)
#%%
def load_IC50(cancer_type, cell_line, IC50_only=True):
    print('TODO: filter out for dosage?')
    df=pd.read_excel(DATA_DIR+'BinChen2017'+os.sep+'SD8.xlsx',\
                     sheet_name=cancer_type, header =0, usecols=['pert_iname', 'pert_id','standard_value', 'standard_units',\
                                 'standard_type', 'cell_line','activity', 'standard_value_median'])
    #filter for cell line
    df = df[df.cell_line==cell_line]
    
    if IC50_only:
        df[df.standard_type=='IC50']
    
    return df.sort_values(by='standard_value_median')

    #filter for 1 nM
#%%   
cancer_type='BRCA'
cell_line='MCF7'
ic50=load_IC50(cancer_type, cell_line)

#%%

ic50.rename( columns = {'pert_iname':'drug'}, inplace=True)
#%%
# calculate spearman coefficient between
#  cell_line cells drug rankings vs IC50s

# merge on common drugs
merged = pd.merge(dr, ic50, on="drug", how="inner", suffixes=("_dr","_ic50"))

merged = merged.dropna(subset=["connectivity_score","standard_value_median"]).copy()
# Following Bin Chen 2017, more negative reversal score 
# (RGES/CS) should be associated with stronger efficacy (lower IC50).
# We report rho for (CS vs IC50) and (CS vs log10(IC50)).
merged["log10_ic50"] = np.log10(merged["standard_value_median"])

rho_linear, p_linear = spearmanr(merged["connectivity_score"], merged["standard_value_median"])
rho_log, p_log = spearmanr(merged["connectivity_score"], merged["log10_ic50"])

print('rho, pval linear:', np.round(rho_linear,2), np.round(p_linear, 2))
print('rho, pval log:', np.round(rho_log, 2), np.round(p_log,2))
