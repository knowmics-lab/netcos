# -*- coding: utf-8 -*-
import os, sys
HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))  # go up 3 levels
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from datetime import datetime
import numpy as np
import pandas as pd
from scipy import stats
import time
import pickle
from loader import load_disease_signature
from connectivity_score import get_common_genes, bin_chen_connectivity
from conf import  BASE_DIR, cell_line, pert_time, DISEASE, MITH_IN_DISEASE, MITH_IN_DRUG, map_name_to_id, MITH_APP, MITH_OUT_DRUG, MITH_OUT_DISEASE,\
CS_IN_DRUG, CS_IN_DISEASE, cell_lines_chembl

print(cell_line, pert_time, DISEASE)
from cs_batch import run_connectivity_score_drugs_batch
from preprocessing_utils import get_drugs_list_from_path, get_chunk_indexes
from conf import DISEASE, CS_OUT, CS_IN_DRUG, CS_IN_DISEASE
print(DISEASE, cell_line, pert_time)
print(CS_IN_DRUG, CS_IN_DISEASE, CS_OUT)

#%%

def load_single_drug_signature(filename, mith=1):
    with open(filename+'.pkl', 'rb') as f:
        data=pickle.load(f)
    return data

def run_connectivity_score_drugs_batch(DISEASE, mith, drugs_list, pert_times, i1, i2, rank_on='magnitude', save_file=False):
    '''
    wrapper function to load disease and drug data (for a  batch of drugs in drugs_list[i1:i2]), 
    and calculate their connectivity score, for given drugs and LINCS drug perturbation times.
    Also calculates Pearson correlation
    coefficient, Spearman correlation coefficient, and cosine similarity, between
    disease and drug.
    Output: a dataframe of connectivity scores and  pvalues, for all combinations of conditions
    '''
    start=time.time()
    
    singature_uom = 'DE_log2_FC' if not mith else 'Perturbation'
    
    # load disease signature:
    disease_signature=load_disease_signature(DISEASE, mith=mith)
    disease_signature=disease_signature[['gene_id', singature_uom, 'adj.p.value']]
    
    
    # Discard disease genes that are 
    # not found in drug genes
    # DO NOT simply select common genes, as this would impact
    # connectrivity calculations
        #load one drug signature (all durg signatures have the same genes in the same order)
    example_drug = drugs_list[0]
    example_drug_filename=CS_IN_DRUG+example_drug
    drug_signature=load_single_drug_signature(example_drug_filename, mith=mith)
    
    disease_common_index, drug_common_index = get_common_genes(disease_signature, drug_signature)  

    
    # filter disease signature with common genes:
    # Note: indexing every time with pandas .iloc is slower than writing
    # and loading binaries for filtered dataframes, for each drug
    # however this is not relevant, as we only need to filter
    # disease data, while we can take all drug data.
    disease_signature=disease_signature.iloc[disease_common_index].reset_index(drop=True)
    print('common disease genes', len(disease_signature))
    

    
    # Initialize dataframe:
    data = []
    
    
    # Calculate connectivity score between disease and drugs:     
    for drug in drugs_list[i1:i2]:
        print(drug)
        
        # load drug signature
        
        drug_signature=load_single_drug_signature(CS_IN_DRUG+drug, mith=mith)
        
        for pert_time in pert_times:
            
            # Calculate connectivity score between drug and disease:
            
            # select columns of interest:
            
            p_val_str= 'adj.p.value' if not pert_time=='6h_24h' else 'p.value'
            columns_of_interest=['gene_id', singature_uom+'_'+pert_time, p_val_str+'_'+pert_time]
            cscore, p_value = bin_chen_connectivity(disease_signature, drug_signature[columns_of_interest], rank_on=rank_on)
            
            # Calculate other correlations (between common genes)
            disease_signature_signature=disease_signature[singature_uom]
            drug_signature_signature=drug_signature[singature_uom+'_'+pert_time].loc[drug_common_index]
            pearson=stats.pearsonr(disease_signature_signature, drug_signature_signature)
            spearman=stats.spearmanr(disease_signature_signature, drug_signature_signature)
            cosine_sim=np.dot(np.array(disease_signature_signature),np.array(drug_signature_signature))/(np.linalg.norm(np.array(disease_signature_signature))*np.linalg.norm(np.array(drug_signature_signature)))
            
            data.append([DISEASE, drug, pert_time, cscore, p_value, pearson[0], pearson[1], spearman[0], spearman[1], cosine_sim]) # other correlations, genes subset
                
    
    connectivity_data= pd.DataFrame(data, columns=["disease","drug","perturbation_time","connectivity_score","cs_p_value",'pearson','pearson_p_value','spearman','spearman_p_value','cos_sim']) #add other correlations, genes subset
    
    # save connectivity scores
    if save_file:
        if not os.path.exists(CS_OUT):
            os.mkdir(CS_OUT)
        connectivity_dataset_filename=CS_OUT+str(i1)+'_'+str(i2)+'_DEG_connectivity_score.tsv' if not mith else  CS_OUT+str(i1)+'_'+str(i2)+'_mith_connectivity_score.tsv'
        connectivity_data.to_csv(connectivity_dataset_filename, sep='\t', index=False)
    
    print('total elapsed time for batch of', len(drugs_list[i1:i2]),' drugs: ', time.time()-start)
    return connectivity_data
#%%
mith=1
drugs_list=get_drugs_list_from_path(CS_IN_DRUG)
print(len(drugs_list))
# set n of cores for parallel execution
n_jobs= 4
if n_jobs>len(drugs_list):
    raise ValueError('n_jobs:',n_jobs,' > len(drugs_list):',len(drugs_list),'! Reduce size of n_jobs')

# get chunk size
chunk_size=int(len(drugs_list)/n_jobs)
last_chunk_size=len(drugs_list)%n_jobs

pert_times=[pert_time]

parallel_indexes=[]
for i in range(n_jobs):
    i1,i2=get_chunk_indexes(i, chunk_size)
    parallel_indexes.append((i1,i2))

print('n jobs', n_jobs)
print('len drug list', len(drugs_list))
print('chunk size',chunk_size)
print('last chunk size', last_chunk_size)
#%%
run_connectivity_score_drugs_batch(DISEASE, mith, drugs_list, pert_times, i1, i2, save_file=True)
