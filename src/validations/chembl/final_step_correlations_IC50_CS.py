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
HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))  # go up 3 levels
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from conf import DISEASE, CS_OUT, DATA_DIR, CS_DIR,\
    pert_time, cell_line, diseases_of


# select 
def load_drug_rankings(path, filename=None, pert_time='all',mith='mith'):
    
    
    if not filename:
        filename=mith+'_connectivity_score.tsv'
    df= pd.read_csv(path+filename, sep='\t', header=0, usecols=[0,1,2,3])
    
    if not pert_time == 'all':
        return dr[dr.perturbation_time==pert_time]
    
    return df

def translate_cl(cell_line):
    if cell_line == 'HT29':
        return 'HT-29'
    return cell_line

def load_IC50(cancer_type, cell_line, IC50_only=True):
    
    cell_line = translate_cl(cell_line)
    df=pd.read_excel(DATA_DIR+'BinChen2017'+os.sep+'SD8.xlsx',\
                     sheet_name=cancer_type, header =0, usecols=['pert_iname', 'pert_id','standard_value', 'standard_units',\
                                 'standard_type', 'cell_line','activity', 'standard_value_median'])
    #filter for cell line
    df = df[df.cell_line.str.upper()==cell_line]
    
    if IC50_only:
        df[df.standard_type=='IC50']
    
    # rename for later merge
    df.rename( columns = {'pert_iname':'drug'}, inplace=True)
    return df.sort_values(by='standard_value_median')

#%%
if __name__=="__main__":
    
    print('running correlations between chembl IC50 and drug rankings for disease:', DISEASE)
    print(DISEASE, CS_OUT, pert_time, cell_line)
    dr=load_drug_rankings(CS_OUT)
    print(dr.shape, 'drugs ')
    #%%
    
        #filter for 1 nM
    #%%   
    ic50=load_IC50(DISEASE, cell_line)
    print(ic50.shape, 'drugs with IC50 value for cell line', cell_line)
    #%%
    
    
    #%% merge on common drugs
    merged = pd.merge(dr, ic50, on="drug", how="inner", suffixes=("_dr","_ic50"))
    
    merged = merged.dropna(subset=["connectivity_score","standard_value_median"]).copy()
    print('merged data on common drugs:', merged.shape)
    
    #%% Optional: check overlap between our overlap of ic50 vs CS score
    # and BinChen's overlap of ic50vs sRGES:
    
    BC_merged = df=pd.read_excel(DATA_DIR+'BinChen2017'+os.sep+'SD5.xlsx',\
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
    
    rho_linear, p_linear = spearmanr(BC_merged["sRGES"], BC_merged["standard_value"])
    rho_log, p_log = spearmanr(BC_merged["sRGES"], BC_merged["log10_ic50"])
    
    print('linear IC50: rho=', np.round(rho_linear,2),' pval =' ,np.round(p_linear, 2))
    print('log IC50: rho=', np.round(rho_log, 2),' pval =', np.round(p_log,2))
    
    #%%
    # da fare dopo pranzo: raccogliere i calcoli di sopra in un dizionario che
    # abbia tutte le voci, e riempirlo per tutte le disease e eprt time
    # renderlo poi un dizionario, e traslare tutto su ntobeook (lo stesso di cehmbl validation)
    # per printare il dataframe per bene. 
    
    # salvare eventualmente in un dizionario a parte,
    # i valori di IC50_drugs, common_drugs, overlap Drugs
    data_of = {}
    for cell_line in ['MCF7','HEPG2','HT29']:
        print('-----------------------',cell_line)
        disease = diseases_of[cell_line]
        ic50 = load_IC50(disease, cell_line)
        BC_merged = df=pd.read_excel(DATA_DIR+'BinChen2017'+os.sep+'SD5.xlsx',\
                                     sheet_name=disease)
        print(ic50.shape, 'drugs with IC50 value for cell line', cell_line)
    
        for pert_time in ['6h', '24h']:
            cs_out=CS_DIR+os.sep+'output'+os.sep+disease+'_2025_'+pert_time+os.sep
            print(cell_line, 'pert time',pert_time)
            dr=load_drug_rankings(cs_out)
    
            merged = pd.merge(dr, ic50, on="drug", how="inner", suffixes=("_dr","_ic50"))
            merged = merged.dropna(subset=["connectivity_score","standard_value_median"]).copy()
            print('merged data on common drugs:', merged.shape)
            
            overlap_BC = set(merged.drug).intersection(set(BC_merged.pert_iname))
            print('overlap between merged data with MCS vs merged data with sRGES from BinCHen2017', len(overlap_BC))
            
            merged["log10_ic50"] = np.log10(merged["standard_value_median"])
    
            rho_linear, p_linear = spearmanr(merged["connectivity_score"], merged["standard_value_median"])
            rho_log, p_log = spearmanr(merged["connectivity_score"], merged["log10_ic50"])
    
            print('linear IC50: rho=', np.round(rho_linear,2),' pval =' ,np.round(p_linear, 2))
            print('log IC50: rho=', np.round(rho_log, 2),' pval =', np.round(p_log,2))
            
            BC_merged["log10_ic50"] = np.log10(BC_merged["standard_value"])
            
            BC_rho_linear, BC_p_linear = spearmanr(BC_merged["sRGES"], BC_merged["standard_value"])
            BC_rho_log, BC_p_log = spearmanr(BC_merged["sRGES"], BC_merged["log10_ic50"])
            
            print('linear IC50: rho=', np.round(rho_linear,2),' pval =' ,np.round(p_linear, 2))
            print('log IC50: rho=', np.round(rho_log, 2),' pval =', np.round(p_log,2))
            
            data_of[disease+'_'+pert_time] = [ic50.shape[0], dr.shape[0],\
                                              merged.shape[0],BC_merged.shape[0],\
                                            len(overlap_BC),rho_linear, p_linear,\
                                                rho_log, p_log, BC_rho_linear,\
                                                BC_p_linear, BC_rho_log, BC_p_log]
                
    
    data = pd.DataFrame(data_of).transpose()
    data.columns= ['n_IC50','n_MCS', 'n_MCSvIC50', \
                                            'n_sRGESvIC50','n_overlapSRGESvMCS',\
                                            'rho','p','rho_log','p_log',\
                                            'rho_sRGES', 'p_sRGES','rho_log_sRGES',\
                                                'p_log_sRGES']
