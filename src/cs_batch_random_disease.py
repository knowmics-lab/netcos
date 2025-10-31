# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 12:08:22 2025

@author: L-F-S

Parallel calculation of conectivity score and other similarity measures
for a list of drugs vs a list of random disease gene rankings
#%% Altra idea, invece di vedere quanto il mnostro drug ranking si distanzi
# da un ranking di disease casuale, quindi confrontando le nostre 3k drug FC oppure
# drug PF contro tipo 10k ranking casuali di geni. Quindi prendere le FC vere
# e confrontarle contro drug rankings, a questo punto potrei fare delle correlazioni
# xke cosi ho un drug ranking sulle drug esistenti vere e proprie vs un disease c
# casuale. Questa sicuramente da mettere in un altro notebook


"""

import os
import sys
import numpy as np
import pandas as pd
import time
from connectivity_score import get_common_genes, random_disease_bin_chen_connectivity
from loader import load_disease_signature, load_single_drug_signature
from conf import DISEASE, CS_OUT
from preprocessing_utils import get_drugs_list, get_chunk_indexes
from joblib import Parallel, delayed


###############################################################################
#             CALCULATE CONNECTIVITY FOR DEG/MITH DATA
###############################################################################
def run_disease_randomized_connectivity_score_drugs_batch(random_disease_list, mith, drugs_list, i1, i2, rank_on='magnitude', save_file=False):
    '''
    wrapper function to load disease and drug data (for a  batch of drugs in drugs_list[i1:i2]), 
    and calculate their connectivity score, for given drugs and LINCS drug perturbation times.
    Also calculates Pearson correlation
    coefficient, Spearman correlation coefficient, and cosine similarity, between
    disease and drug.
    inputs:
        random_disease: str: rndom disease name
    Output: a dataframe of connectivity scores and  pvalues, for all combinations of conditions
    '''
    start=time.time()
    
    singature_uom = 'DE_log2_FC' if not mith else 'Perturbation'
    
    # Initialize dataframe:
    data = []
    
    
    # Calculate connectivity score between disease and drugs: 
    for random_disease in random_disease_list[i1:i2]:
        random_disease_name='RD'+str(random_disease)
        for drug in drugs_list:
            print(drug)
            
            # load drug signature
            
            drug_signature=load_single_drug_signature(drug, mith=mith)
            
                
            # Calculate connectivity score between drug and disease:
            pert_time='6h_24h'
            columns_of_interest=['gene_id', singature_uom+'_'+pert_time]
            cscore = random_disease_bin_chen_connectivity(drug_signature[columns_of_interest], rank_on=rank_on)
            data.append([random_disease_name, drug, pert_time, cscore]) # other correlations, genes subset
                    
    connectivity_data= pd.DataFrame(data, columns=["disease","drug","perturbation_time","connectivity_score"]) #add other correlations, genes subset
    
    # save connectivity scores
    if save_file:
        if not os.path.exists(CS_OUT):
            os.mkdir(CS_OUT)
        connectivity_dataset_filename=CS_OUT+str(i1)+'_'+str(i2)+'RD_DEG_connectivity_score.tsv' if not mith else  CS_OUT+str(i1)+'_'+str(i2)+'RD_mith_connectivity_score.tsv'
        connectivity_data.to_csv(connectivity_dataset_filename, sep='\t', index=False)
    
    print('total elapsed time for batch of', len(random_disease_list[i1:i2]),' random diseases: ', time.time()-start)
    return connectivity_data

#%% Parallel run
# get drugs list
# set perturbation times
# Set calculation for MITHrIL data
mith=1
drugs_list=get_drugs_list(mith)
n_random_diseases=10

# set n of cores for parallel execution
n_jobs= int(sys.argv[1]) #16
if n_jobs>n_random_diseases:
    raise ValueError('n_jobs:',n_jobs,' > n_random_diseases:',n_random_diseases,'! Reduce size of n_jobs')

# get chunk size
chunk_size=int(n_random_diseases/n_jobs)
last_chunk_size=n_random_diseases%n_jobs

parallel_indexes=[]
for i in range(n_random_diseases):
    i1,i2=get_chunk_indexes(i, chunk_size)
    parallel_indexes.append((i1,i2))

print('n jobs', n_jobs)
print('len drug list', len(drugs_list))
print('n_random_diseases', n_random_diseases)
print('chunk size',chunk_size)
print('last chunk size', last_chunk_size)
print('mith tag:', mith, '\ndisease:', DISEASE)
results = Parallel(n_jobs=n_jobs)(delayed(run_disease_randomized_connectivity_score_drugs_batch)\
                             (np.arange(n_random_diseases), mith, drugs_list,  i1, i2)\
                        for i1,i2 in parallel_indexes)

if last_chunk_size>0:
    last_batch=run_disease_randomized_connectivity_score_drugs_batch(np.arange(n_random_diseases), mith, drugs_list, i2,n_random_diseases)
    results.append(last_batch)

# build dataframe    
cs_df=pd.concat(results)

# write results
if not os.path.exists(CS_OUT):
            os.mkdir(CS_OUT)
connectivity_dataset_filename=CS_OUT+'RD_DEG_connectivity_score.tsv' if not mith else  CS_OUT+'RD_mith_connectivity_score.tsv'
#%%
cs_df.to_csv(connectivity_dataset_filename, sep='\t', index=False)
