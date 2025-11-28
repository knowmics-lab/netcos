# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 09:08:11 2025

@author: los4
"""
import pandas as pd
import os
import pyreadr
current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(os.path.dirname(os.path.dirname(current_file_path)))
import sys
sys.path.insert(0, parent_directory)
from conf import TSR_OUT_DRUG, cell_line, pert_time,\
    MITH_IN_DRUG
DRUG_DIR = TSR_OUT_DRUG+'LINCS'+os.sep+cell_line+'_'+pert_time+'_drug_wise'+os.sep
df_list=[]
for i, drug_file in enumerate(os.listdir(DRUG_DIR)):
    drug = drug_file.rstrip('.Rds')
    result = pyreadr.read_r(DRUG_DIR+drug_file)
    print(drug)
    df = result[None] 
    df.set_index('gene_id', inplace=True)
    df.rename(columns = {'DE_log2_FC':drug}, inplace = True)
    
    df_list.append(df[drug])

matrix=pd.concat(df_list, axis=1)
print(matrix.shape)
print(matrix.columns)
print(matrix.head(3))

# Save MITHrIL gene id based input file:
mith_input_filename = 'LINCS_'+cell_line+'_'+pert_time+'.mi'
matrix.to_csv(MITH_IN_DRUG+mith_input_filename, sep='\t', index=True,  index_label=False)
print('MITHrIL drug input saved in',MITH_IN_DRUG,'filename:',mith_input_filename)