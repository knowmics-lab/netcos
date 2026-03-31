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

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from conf import DISEASE, CS_OUT, DATA_DIR, CS_DIR,\
    cell_line, diseases_of, LOGS_DIR,\
        cs_filename, disease_run_name, cell_line_run_name,\
    cs_on_LM, cs_mith, selected_cs_run_id, \
    cs_log_filename, lincs_metadata_path, chembl_val_log_filename
from logger import append_run_metadata


def resolve_cs_run_id(
    cs_runs_tsv,
    disease_run_id,
    drug_run_id,
    cs_on_LM,
    mith,
    selected_cs_run_id=None,
):
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

def add_pert_id_to_cs( lincs_metadata_path,cs_df,
                      cs_id_col='LINCS_id',
                      metadata_id_col='id',
                      metadata_pert_col='pert_id'):
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
# select 
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

def load_IC50(ic50_file, cancer_type, cell_line, IC50_only=True):
    
    cell_line = translate_cl(cell_line)
    df=pd.read_excel(ic50_file,\
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
    
    dr=load_drug_rankings(CS_OUT, filename = cs_drug_file)
    print(dr.shape, 'drugs ')
    #%%
    
        #filter for 1 nM
    #%%
    ic50_file = DATA_DIR/'BinChen2017'/'SD8.xlsx'
    ic50=load_IC50(ic50_file, DISEASE, cell_line)
    print(ic50.shape, 'drugs with IC50 value for cell line', cell_line)
    
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
    
    print('linear IC50 SD5: rho=', np.round(rho_linear,2),' pval =' ,np.round(p_linear, 2))
    print('log IC50 SD5: rho=', np.round(rho_log, 2),' pval =', np.round(p_log,2))
    
    #%%
    


    
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
        "n_cs_rows": len(dr),
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

        # output
        # "output_plot_file": str(output_plot_file)
    }
    

    append_run_metadata(chembl_val_log_filename, run_metadata_data)
    
    
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
