#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 09:31:28 2025

@ author: L-F-S

utility tests
"""

import os

import sys
import numpy as np
import pandas as pd
import pickle
import time
from loader import load_single_drug_signature, load_drug_signatures
from conf import BASE_DIR, DISEASE,MITH_IN_DRUG, MITH_OUT_DRUG,\
 TSR_OUT_DRUG, TSR_OUT_DISEASE, TSR_OUT_CSCORE, alias_2geneid
 #%% test pickled
# n=0
# for df_file in os.listdir(TSR_OUT_DRUG+'/LINCS/metanalysis_drug_wise_filtered/'):
#     if df_file.endswith('.pkl'):
#         drug=df_file.split('_')[0]
#         n+=1
#         if n%100==0:
#             print(n)
#         data_pkl=load_single_drug_signature(drug, mith=False, pkl=True)
#         data=load_single_drug_signature(drug, mith=False, pkl=False)
        
#         same=np.unique(data==data_pkl)[0]
#         if not same:
#             print('WARNING test failed for drug', drug)
        
# n=0  
# for df_file in os.listdir(TSR_OUT_DRUG+'/LINCS/metanalysis_mith3_drug_wise/'):
#     if df_file.endswith('.csv'):
#         drug=df_file.split('_')[0]
#         n+=1
#         if n%100==0:
#             print(n)
#         data_pkl=load_single_drug_signature(drug, mith=False, pkl=True)
#         data=load_single_drug_signature(drug, mith=False, pkl=False)
        
#         same=np.unique(data==data_pkl)[0]
        
#         if not same:
#             print('WARNING test failed for drug', drug)

# print('end')

#%% test load speed

start=time.time()
load_drug_signatures(pkl=True)
print(time.time()-start, 'to load pickled')
# 82.10381436347961 to load pickled for 3222 drugs

start=time.time()
load_drug_signatures(pkl=False)
print(time.time()-start, 'to load unpickled')  
# 286.05032181739807 to load unpickled for 3222 drugs

#%%% TESTS for connectivity score speed optimization
# def montecarlo_connectivity(s_up, s_down, r, n_iterations=1000):
#     '''
#     Randomly sample RGES
#     '''
#     random_RGES_list=[]
#     disease_genes=range(s_up+s_down) # the index are the sorted fake disease genes
#     for i in range(n_iterations):
        
#         # Generate random indexes for up and downregulated genes
#         # random_up_and_down_indexes=random.sample(range(r), s_up+s_down) using python's random (slower)
#         # fastest implementation using numpy:
#         random_up_and_down_indexes=np.random.choice(r, s_up + s_down, replace=False)
         
         
#         # Create random index dictionary:
            
#         # random_V=pd.Series(index=disease_genes, data=random_up_and_down_indexes).to_dict()
#         # random_V={i:random_up_and_down_indexes[i] for i in disease_genes}
        
#         # Fastest implementation:
#         random_V=dict(zip(disease_genes, random_up_and_down_indexes))
        
#         # TEST to compare they are indeed all the same
#         # print('comparing')
#         # for i1, i2 in random_V.items():
#         #     if not random_V2[i1]==i2:
#         #         print(i1,i2)
#         # print('done')
        
#         # Compute random KS stats
#         random_a_up, random_b_up, _, _ = compute_Kolmogorov_Smirnof_stat(random_up_and_down_indexes[:s_up], random_V)
#         random_a_down , random_b_down, _, _ = compute_Kolmogorov_Smirnof_stat(random_up_and_down_indexes[-s_down:], random_V)
        
#         # Calculate random RGES
#         random_RGES_list.append(calculate_RGES(random_a_up, random_a_down, random_b_up, random_b_down))
     
#     return random_RGES_list

# def rank_genes(drug_signature, disease_signature, drug_col_name, disease_col_name):
#     '''
#     output:
#         dictionary where keys are sorted disease gene indexes 
#         and values are indexes of sorted drug genes
        
#     '''
    
#     sorted_drug_genes=drug_signature.sort_values(by=drug_col_name).reset_index(drop=True) 
#     sorted_disease_genes=disease_signature.sort_values(by=disease_col_name).reset_index(drop=True)
#     merged_data_on_disease_indexes = pd.merge(sorted_disease_genes, sorted_drug_genes.reset_index(), on='gene_id')
#     V=merged_data_on_disease_indexes['index'].to_dict() # create a dictionary of sorted disease index : sorted drug index
    
#     # TEST: manually iterate over disease indexes and build dictionary, and compare they are equal (they are)
#     # V2={i:sorted_drug_genes[sorted_drug_genes.gene_id==sorted_disease_genes.iloc[i].gene_id].index[0] for i in sorted_disease_genes.index}
#     # print('comparing')
#     # for i1, i2 in V.items():
#     #     if not V2[i1]==i2:
#     #         print(i1,i2)
#     # print('done')
#     return V

def bin_chen_connectivity(disease_signature, drug_signature):
    '''calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017
   Inputs:
          disease_signature: pd.DataFrame() columns: gene_id, signature data (DE_log2_FC or mith perturbation), p value
          drug_signature: pd.DataFrame() columns: gene_id, signature data (DE_log2_FC or mith perturbation), p value
    '''
    disease_p_val_col_name=disease_signature.columns[2]
    disease_signature_col_name=disease_signature.columns[1]
    
    drug_p_val_col_name=drug_signature.columns[2]
    drug_signature_col_name=drug_signature.columns[1]

    # # Get lists of up (down) regulated genes: older with 2vs
    disease_signature_up = disease_signature[disease_signature[disease_signature_col_name]>0]
    disease_signature_down = disease_signature[disease_signature[disease_signature_col_name]<0]
        
    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, disease_signature_up  , drug_p_val_col_name, disease_p_val_col_name)
    V_down = rank_genes(drug_signature, disease_signature_down, drug_p_val_col_name, disease_p_val_col_name)
#    V = rank_genes(drug_signature, disease_signature, drug_p_val_col_name, disease_p_val_col_name)

    # Compute KS statistic for up and down regulated genes:
    # # Get lists of up (down) regulated genes:
    # disease_up_i = disease_signature[disease_signature[disease_signature_col_name]>0].index
    # disease_down_i = disease_signature[disease_signature[disease_signature_col_name]<0].index
    # a_up, b_up, s_up, r = compute_Kolmogorov_Smirnof_stat(disease_signature['gene_id'].loc[disease_up_i], V, drug_signature['gene_id'])
    # a_down , b_down, s_down, r = compute_Kolmogorov_Smirnof_stat(disease_signature['gene_id'].loc[disease_down_i], V, drug_signature['gene_id'])
    a_up, b_up, s_up, r = compute_Kolmogorov_Smirnof_stat(disease_signature_up['gene_id'], V_up, drug_signature['gene_id'])
    a_down , b_down, s_down, r = compute_Kolmogorov_Smirnof_stat(disease_signature_down['gene_id'], V_down, drug_signature['gene_id'])


    # Compute RGES:
    measured_RGES = calculate_RGES(a_up, a_down, b_up, b_down)
    
    # Calculate two tailed p-value (no assumption on the direction of RGES 
    # between disease and drug) for measured RGES, using random sampling:
    n_iterations=1000 
    random_RGES_list=montecarlo_connectivity(s_up, s_down, r, n_iterations)
    p_value=np.sum(np.abs(np.array(random_RGES_list))>np.abs(measured_RGES))/n_iterations

    return measured_RGES, p_value

#%% test vectorized version of connectivity score has the same results
def rank_genes(drug_signature, disease_signature, drug_col_name, disease_col_name):
    '''
    output:
        dictionary where keys are sorted disease gene indexes 
        and values are indexes of sorted drug genes
        
    '''
    #TODO account for 1 based indexing if needed
    
    sorted_drug_genes=drug_signature.sort_values(by=drug_col_name).reset_index(drop=True) 
    sorted_disease_genes=disease_signature.sort_values(by=disease_col_name).reset_index(drop=True)
    merged_data_on_disease_indexes = pd.merge(sorted_disease_genes, sorted_drug_genes.reset_index(), on='gene_id')
    V=merged_data_on_disease_indexes['index'].to_dict() # create a dictionary of sorted disease index : sorted drug index
    
    return V

def compute_Kolmogorov_Smirnof_stat(disease_disregulated_genes, V, drug_genes):
    '''From Lamb et al., 2006 supplementary'''
    
    # number of disregulated genes:
    s=len(disease_disregulated_genes) # TODO called num_tags_up(down) in Catalano,
                                      # there it is the length of the disease 
                                      # genes+the drug vector genes.
    
    # total number of genes in reference drug expression data:
    r=len(drug_genes) #num_genes in catalano
    
    # Compute statistic (accounting for python's 0 based indexing)
    a=max([(j+1)/s-V[j]/r for j in range(s)]) 
    b=max([V[j]/r-j/s for j in range(s)])
    
    return a, b, s, r

def vectorized_rank_genes(drug_signature, disease_signature, drug_col_name, disease_col_name):
    '''
    output:
        NumPy array where the indices are sorted disease gene indexes 
        and the values are indexes of corresponding sorted drug genes
    '''
    
    sorted_drug_genes=drug_signature.sort_values(by=drug_col_name).reset_index(drop=True)
    sorted_disease_genes=disease_signature.sort_values(by=disease_col_name)
    merged_data_on_disease_indexes = pd.merge(sorted_disease_genes, sorted_drug_genes.reset_index(), on='gene_id')
    
    V=merged_data_on_disease_indexes['index'].to_numpy() 
    return V

def vectorized_compute_KS(disease_disregulated_genes, V, drug_genes):
    '''From Lamb et al., 2006 supplementary'''
    
    # number of disregulated genes:
    s=len(disease_disregulated_genes) 

    # total number of genes in reference drug expression data:
    r=len(drug_genes) #num_genes in catalano

    
    V_over_r = V / r
    
    # Calculate a and b using vectorized operations
    a = np.max(np.arange(1,s+1) / s - V_over_r)
    b = np.max(V_over_r - (np.arange(s) / s))
    
    return a, b, s, r

def vectorized_montecarlo_connectivity(s_up, s_down, r, n_iterations=1000):
    '''
    Randomly sample RGES
    '''
    random_RGES_list=[]
    disease_genes=range(s_up+s_down) # the index are the sorted fake disease genes
    drug_genes=np.arange(r)
    
    for i in range(n_iterations):
        
        # Generate random indexes for up and downregulated genes
        # fastest implementation using numpy:
        random_up_and_down_indexes=np.random.choice(r, s_up + s_down, replace=False)
        # slower:
        # random_up_and_down_indexes=random.sample(range(r), s_up + s_down)
         
        # Create random index dictionary:
        # Fastest implementation:
        random_V_up=random_up_and_down_indexes[:s_up]
        random_V_down=random_up_and_down_indexes[-s_down:]
        
        random_V_up_non=dict(zip(range(s_up), random_up_and_down_indexes[:s_up]))
        random_V_down_non=dict(zip(range(s_down), random_up_and_down_indexes[-s_down:]))
        
        if not list(random_V_up) == list(random_V_up_non.values()):
            raise ValueError('wrong lists are different')

        # Compute random KS stats:
        random_a_up_non, random_b_up_non, _, _ = compute_Kolmogorov_Smirnof_stat(random_up_and_down_indexes[:s_up], random_V_up_non, drug_genes)
        random_a_down_non , random_b_down_non, _, _ = compute_Kolmogorov_Smirnof_stat(random_up_and_down_indexes[-s_down:], random_V_down_non, drug_genes)

        
        # Compute random KS stats:
        random_a_up, random_b_up, _, _ = vectorized_compute_KS(random_up_and_down_indexes[:s_up], random_V_up, drug_genes)
        random_a_down , random_b_down, _, _ = vectorized_compute_KS(random_up_and_down_indexes[-s_down:], random_V_down, drug_genes)
        
        if not float(random_a_up) == float(random_a_up_non):
            raise ValueError('wrong', float(random_a_up), float(random_a_up_non))
        if not float(random_b_up) == float(random_b_up_non):
            raise ValueError('wrong', float(random_b_up), float(random_b_up_non))
        
        # Calculate random RGES:
        random_RGES_list.append(calculate_RGES(random_a_up, random_a_down, random_b_up, random_b_down))
     
    return random_RGES_list

def vectorized_bin_chen_connectivity(disease_signature, drug_signature):
    '''calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017
   Inputs:
          disease_signature: pd.DataFrame() columns: gene_id, signature data (DE_log2_FC or mith perturbation), p value
          drug_signature: pd.DataFrame() columns: gene_id, signature data (DE_log2_FC or mith perturbation), p value
    '''
    disease_p_val_col_name=disease_signature.columns[2]
    disease_signature_col_name=disease_signature.columns[1]
    
    drug_p_val_col_name=drug_signature.columns[2]
    drug_signature_col_name=drug_signature.columns[1]
    
    # # Get lists of up (down) regulated genes: older with 2vs
    disease_signature_up = disease_signature[disease_signature[disease_signature_col_name]>0]
    disease_signature_down = disease_signature[disease_signature[disease_signature_col_name]<0]
    
    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, disease_signature_up  , drug_p_val_col_name, disease_p_val_col_name)
    V_down = rank_genes(drug_signature, disease_signature_down, drug_p_val_col_name, disease_p_val_col_name)
    
    V_up_vector = vectorized_rank_genes(drug_signature, disease_signature_up  , drug_p_val_col_name, disease_p_val_col_name)
    V_down_vector = vectorized_rank_genes(drug_signature, disease_signature_down, drug_p_val_col_name, disease_p_val_col_name)
   
    if not list(V_up_vector) == list(V_up.values()):
        raise ValueError('wrong lists are different', V_up_vector[:10], list(V_up.values())[:10] )
    
    # Compute KS statistic for up and down regulated genes:
    a_up_unv, b_up, s_up, r = compute_Kolmogorov_Smirnof_stat(disease_signature_up['gene_id'], V_up, drug_signature['gene_id'])
    a_down , b_down_unv, s_down, r = compute_Kolmogorov_Smirnof_stat(disease_signature_down['gene_id'], V_down, drug_signature['gene_id'])

    a_up, b_up, s_up, r = vectorized_compute_KS(disease_signature_up['gene_id'], V_up_vector, drug_signature['gene_id'])
    a_down , b_down, s_down, r = vectorized_compute_KS(disease_signature_down['gene_id'], V_down_vector, drug_signature['gene_id'])
    
    if not float(a_up) == float(a_up_unv):
        raise ValueError('wrong', float(a_up), float(a_up_unv))
    if not float(b_down) == float(b_down_unv):
        raise ValueError('wrong', float(b_down), float(b_down_unv))
        
    
    # Compute RGES:
    measured_RGES = calculate_RGES(a_up, a_down, b_up, b_down)
    
    
    
    # Calculate two tailed p-value (no assumption on the direction of RGES 
    # between disease and drug) for measured RGES, using random sampling:
    n_iterations=1000 
    start=time.time()
    random_RGES_list=vectorized_montecarlo_connectivity(s_up, s_down, r, n_iterations)
    print('vecotirzed',time.time()-start)
    p_value=np.sum(np.abs(np.array(random_RGES_list))>np.abs(measured_RGES))/n_iterations

    return measured_RGES, p_value