#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:41:12 2024

@ author: L-F-S
"""

#%%
from conf import CS_OUT
import pandas as pd
import subprocess
from preprocessing_utils import get_drugs_list
import os

print('aggiungi cosa per levare duplicati')
quit()
#%% Get all drugs
drugs_list=get_drugs_list()

# set mithril flag (0:DEG data, 1:mithril data)
mith=0

# set n of cores
cores=16

# set ranking metric
rank_on='magnitude'

#%%

# get chunk size
chunk_size=int(len(drugs_list)/cores)
last_chunk_size=len(drugs_list)%cores

def get_chunk_indexes(i,chunk_size):
    i1=chunk_size*i
    i2=chunk_size*i+chunk_size
    return i1, i2

#%% merge files
print('merging drug files')
dfs=[]
for i in range(cores):
    i1,i2=get_chunk_indexes(i, chunk_size)
    
    # RUN CS AND SAVE DATA TO FILE
    connectivity_dataset_filename=CS_OUT+str(i1)+'_'+str(i2)+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(i1)+'_'+str(i2)+'_mith_connectivity_score.tsv'
    # LOAD DATA FROM SAVED FILE INTO DATAFRAME
    dfs.append(pd.read_csv(connectivity_dataset_filename, sep='\t'))

# RUN CS AND SAVE DATA TO FILE
connectivity_dataset_filename=CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_mith_connectivity_score.tsv'
# LOAD DATA FROM SAVED FILE INTO DATAFRAME
dfs.append(pd.read_csv(connectivity_dataset_filename, sep='\t'))




connectivity_dataset_filename=CS_OUT+'DEG_connectivity_score_'+rank_on+'.tsv' if mith==0 else CS_OUT+'mith_connectivity_score_'+rank_on+'.tsv'

connectivity_score_df=pd.concat(dfs)

# drop posible duplicated rows for indexing overlap:
    # DA TESTARE
# connectivity_score_df.drop_duplicates(subset=['drug','perturbation_time'], inplace=True)

# #Write data frame
connectivity_score_df.to_csv(connectivity_dataset_filename, sep='\t', index=False)

#%% delete single files
# for i in range(cores):
#     i1,i2=get_chunk_indexes(i, chunk_size)
#     connectivity_dataset_filename=CS_OUT+str(i1)+'_'+str(i2)+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(i1)+'_'+str(i2)+'_mith_connectivity_score.tsv'
#     os.remove(connectivity_dataset_filename)

# connectivity_dataset_filename=CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_mith_connectivity_score.tsv'
# os.remove(connectivity_dataset_filename)

#%%
