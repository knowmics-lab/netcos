#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:10:47 2024

@ author: L-F-S

translate MITHrIL2 into MITHrIL3 input for drug data,
by appending multiple drug columns for all genes
and converting gene names to gene ids
"""
import os
import pandas as pd
from conf import MITH_IN_DRUG, alias_2geneid

df_list=[]
for mi2 in os.listdir(MITH_IN_DRUG):
    drug_condition_name=mi2.rstrip('.mi')
    if not drug_condition_name.startswith('LINCS'):
        print(drug_condition_name)
        df=pd.read_csv(MITH_IN_DRUG+mi2,sep='\t',header=None,index_col=0, names=[drug_condition_name])
        print(df.shape)
        df_list.append(df)

matrix=pd.concat(df_list, axis=1)
print(matrix.shape)
print(matrix.columns)
print(matrix.head(3))

#%%
def map_name_to_id(genename):
    if genename in alias_2geneid.keys():
        return alias_2geneid[genename]
    return genename
matrix.index=matrix.reset_index()['index'].apply(lambda x : map_name_to_id(x))
#%%
LINCS_metanalysis_filename_id='LINCS_metanalysis.mi'
matrix.to_csv(MITH_IN_DRUG+LINCS_metanalysis_filename_id,sep='\t', index=True,  index_label=False) # will save the index but not the index name
print('MITHrIL drug input saved in',MITH_IN_DRUG,'filename:',LINCS_metanalysis_filename_id)
