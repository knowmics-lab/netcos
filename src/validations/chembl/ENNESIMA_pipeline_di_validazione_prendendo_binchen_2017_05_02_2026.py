# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 16:09:56 2026

@author: los4
"""
import os

# this works with conf file created on 18/11 related to chembl validation:
from conf import BASE_DIR, TSR_OUT_DRUG, cell_line, pert_time, DISEASE, MITH_IN_DISEASE, MITH_IN_DRUG, map_name_to_id, MITH_APP, MITH_OUT_DRUG, MITH_OUT_DISEASE,\
CS_IN_DRUG, CS_IN_DISEASE, cell_lines_chembl, diseases_of, DATA_DIR, CS_DIR

from validations.chembl.chembl_loader import load_raw_chembl
import numpy as np
from preprocessing_utils import get_drugs_list
import pandas as pd
import pyreadr
import subprocess
from step3_mith_out_to_cs_in_disease import mith_to_cs_in,write_cs_input
from step6_drug_run_MITHrIL import run_mithril_batch


print(cell_line, pert_time, DISEASE)
from conf import DATA_DIR

DATA_DIR
tcga = pd.read_excel(DATA_DIR+'BinChen2017'+os.sep+'SD2.xlsx', sheet_name = DISEASE+'_sig.csv')
print(tcga.shape, 'all genes')
print(len(tcga[tcga.landmark==1]), 'landmark genes')
#%%
# These three cancer types have
#five (BRCA), two (LIHC) and 12 (COAD) cancer cell lines with
#the same cell lineage (Fig. 1c and Supplementary Data 1) in LINCS data
LINCS_cell_lines_of_disease_df = pd.read_excel(DATA_DIR+'BinChen2017'+os.sep+'SD1.xlsx', sheet_name = DISEASE+'_cell_lines.csv')
LINCS_cell_lines_of_disease=LINCS_cell_lines_of_disease_df.LINCS.dropna()
print(len(LINCS_cell_lines_of_disease))
#%%
