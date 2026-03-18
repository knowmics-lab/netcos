# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 16:09:56 2026

@author: los4
"""
import os
from pathlib import Path

# this works with conf file created on 18/11 related to chembl validation:
from conf import BASE_DIR, TSR_OUT_DRUG, cell_line, pert_time, DISEASE, MITH_IN_DISEASE, MITH_IN_DRUG, map_name_to_id, MITH_APP, MITH_OUT_DRUG, MITH_OUT_DISEASE,\
CS_IN_DRUG, CS_IN_DISEASE, cell_lines_chembl, diseases_of, DATA_DIR, CS_DIR,\
    landmark

from validations.chembl.chembl_loader import load_raw_chembl
import numpy as np
import pandas as pd
import pyreadr

from step6_drug_run_MITHrIL import run_mithril_batch


BC_DATA = DATA_DIR / 'BinChen2017'
LINCS_BC_DATA = BC_DATA / Path('data/data/raw/lincs')

print(' Load LINCS signatures')
result = pyreadr.read_r(LINCS_BC_DATA / 'lincs_signatures_cmpd_landmark.Rdata') # also works for Rds
print(result.keys()) # 'lincs_signatures'

# LINCS signatures between treated and untreated samples, 
# for 978 landmark genes, 
# for 6511 is_gold drugs (at 6h and 24h timepoints)
# rows: landmark gene_id cols: LINCS perturbagen ids
LINCS_FC_bc = result['lincs_signatures']
del result
print(LINCS_FC_bc.shape)
#%%
for cell_line, DISEASE in diseases_of.items():
    print(cell_line, DISEASE, 'landmark?', landmark)
    
    print('load disease data')
    tcga = pd.read_excel(BC_DATA / 'SD2.xlsx', sheet_name = DISEASE+'_sig.csv')
    print(tcga.shape, 'all genes')
    print(len(tcga[tcga.landmark==landmark]), ' genes with LINCS landmark = ' , landmark)
    lm_tcga_fc = tcga[['id','log2FoldChange']][tcga.landmark==landmark]
    lm_tcga_fc['id'] =lm_tcga_fc['id'].apply(lambda x : x.split('|')[1])
    print('convert to mithril input')
    LM_flag = ''
    if landmark:
        LM_flag = '_LM'
    mith_disease_in_filename = DISEASE+LM_flag+'_signature_gene_id.mi'
    lm_tcga_fc.to_csv(MITH_IN_DISEASE/mith_disease_in_filename, sep = '\t', index=False, header=False)
# # These three cancer types have
# # 5 (BRCA), 2 (LIHC) and 12 (COAD) cancer cell lines with
# # the same cell lineage (Fig. 1c and Supplementary Data 1) in LINCS data
# LINCS_cell_lines_of_disease_df = pd.read_excel(DATA_DIR/'BinChen2017'/'SD1.xlsx', sheet_name = DISEASE+'_cell_lines.csv')
# LINCS_cell_lines_of_disease=LINCS_cell_lines_of_disease_df.LINCS.dropna()
# print(len(LINCS_cell_lines_of_disease),LINCS_cell_lines_of_disease)
    print('load ',cell_line,' Binchen2017 FC data')
    
    # gene_id gene_symbol 
    landmark_genes = pd.read_csv(LINCS_BC_DATA / 'lincs_landmark.csv') 
    print(landmark_genes.shape)
    # id for drug at given per time. 66511 drugs (is_gold = 1)
    # relevant columns: id  pert_iname pert_id pert_time cell_id 
    print('load perturbagen metadata')
    drug_md_all = pd.read_csv(LINCS_BC_DATA / 'lincs_sig_info_new.csv') 
    print(drug_md_all.shape)
    print(len(drug_md_all.pert_id.unique()), 'unique compounds')
    
    # filter by cell line
    drug_md = drug_md_all[drug_md_all['cell_id']==cell_line]
    print(drug_md.shape)

    print('filter LINCS signatures by perturbagens appearing in selected cell line')
    
    LINCS_FC_bc_filtered = LINCS_FC_bc[ drug_md.id.astype(str)]
    print(LINCS_FC_bc_filtered.shape)
    
    # 2418 perturbagens for HEPG2

    print('conver LINCS as Mithril input. Remember to manually remove first tab from header!')
    mith_in_lincs_filename = 'LINCS_LM_'+cell_line+'.mi'
    LINCS_FC_bc_filtered.to_csv(MITH_IN_DRUG / mith_in_lincs_filename, header = True , sep = '\t', index = True)
