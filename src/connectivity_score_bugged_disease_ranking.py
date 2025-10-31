#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 20:37:14 2025

@ author: L-F-S
calculate RGES connectivity score between disease signature and drug sinatures.

"""
import os

import sys
import numpy as np
import pandas as pd
import time
from scipy import stats
from conf import DISEASE
from loader import load_disease_signature, load_single_drug_signature

# FUNCTIONS
def get_common_genes(disease_signature, drug_signature):
    '''
    Identifies common genes between drug and disease signature. Returns indexes
    of common genes in both signatures.
    Input:
        - disease_signature: pd.DataFrame with columns ['gene_id', ...]
        - drug_signature: pd.DataFrame with columns ['gene_id', ...]
    Output:
        - Tuple of DataFrame.index of common genes in disease and drug signatures
    
    '''

    missing_genes=set.difference(set(disease_signature.gene_id),set(drug_signature.gene_id))
    # print('missing genes from mithril', len(missing_genes))
    addional_drug_items=set.difference(set(drug_signature.gene_id),set(disease_signature.gene_id))
    # print('addional_mith_items', len(addional_drug_items))
    
    
    common_genes=set.intersection(set(disease_signature.gene_id),set(drug_signature.gene_id))
    # print('genes in common', len(common_genes))
    if not len(common_genes)+len(addional_drug_items)==drug_signature.shape[0]:
        raise ValueError('common genes+ additional mithril items != len(mithril_gene)')
    if not len(common_genes)+len(missing_genes)==disease_signature.shape[0]:
        raise ValueError('riprova: common genes+ missing mithril items != len(DEG_genes)')

    # Get subset of signatures with common genes
    drug_common_genes_signatures=drug_signature[drug_signature['gene_id'].isin(common_genes)].sort_values(by='gene_id')
    disease_common_genes_signatures=disease_signature[disease_signature.gene_id.isin(common_genes)].sort_values(by='gene_id')
    
    if not np.unique(disease_signature.loc[disease_common_genes_signatures.index].reset_index()['gene_id']==drug_signature.loc[drug_common_genes_signatures.index].reset_index()['gene_id']):
        raise ValueError('common gene indexes actually different between drug and disease!')
    return disease_common_genes_signatures.index, drug_common_genes_signatures.index
 

def rank_genes(drug_signature, disease_signature, drug_col_name, disease_col_name):
    '''
    Compute the array V of sorted drug gene indexes, indexed by sorted disease
    indexes.
   Input:
       - drug_signature: pd.DataFrame with columns ['gene_id', drug_col_name, ...]
       - disease_signature: pd.DataFrame with columns ['gene_id', disease_col_name, ...]
       - drug_col_name: str, column name to sort drug_signature by
       - disease_col_name: str, column name to sort disease_signature by
    Output:
       NumPy array where the indices are sorted disease gene indexes 
       and the values are indexes of corresponding sorted drug genes.
   '''
    
    # Sort genes in both datasets:
    sorted_drug_genes=drug_signature.sort_values(by=drug_col_name).reset_index(drop=True)
    sorted_disease_genes=disease_signature.sort_values(by=disease_col_name)
    
    # Merge dataset and keep disease index sorting
    merged_data_on_disease_indexes = pd.merge(sorted_disease_genes, sorted_drug_genes.reset_index(), on='gene_id')
    
    # Extract drug index of sorted genes
    V=merged_data_on_disease_indexes['index'].to_numpy() 
    return V


def compute_KS(disease_disregulated_genes, V, drug_genes):
    '''
    Compute Kolmogorov-Smirnov (KS) statistic.
    From Lamb et al., 2006 supplementary
    Input:
       - disease_disregulated_genes: list or array of gene IDs
       - V: NumPy array of sorted drug gene indices
       - drug_genes: list or array of gene IDs in the drug signature
    Output:
       - a: float, maximum positive difference
       - b: float, maximum negative difference
       - s: int, number of disregulated genes
       - r: int, total number of genes in reference drug expression data
    '''
    
    # number of disregulated genes:
    s=len(disease_disregulated_genes) 

    # total number of genes in reference drug expression data:
    r=len(drug_genes)
    
    V_over_r = V/r
    a = np.max(np.arange(1, s+1)/s - V_over_r)
    b = np.max(V_over_r - (np.arange(s)/s))
    
    return a, b, s, r


def compute_KS_evil_twin(disease_disregulated_genes, V, drug_genes):
    '''
    Compute Kolmogorov-Smirnov (KS) statistic.
    From Lamb et al., 2006 supplementary
    Input:
       - disease_disregulated_genes: list or array of gene IDs
       - V: NumPy array of sorted drug gene indices
       - drug_genes: list or array of gene IDs in the drug signature
    Output:
       - a: float, maximum positive difference
       - b: float, maximum negative difference
       - s: int, number of disregulated genes
       - r: int, total number of genes in reference drug expression data
    '''
    
    # number of disregulated genes:
    s=len(disease_disregulated_genes) 

    # total number of genes in reference drug expression data:
    r=len(drug_genes)
    
    V_over_r = V/r
    a = np.max(np.abs(np.arange(1, s+1)/s - V_over_r))
    b = np.max(np.abs(V_over_r - (np.arange(s)/s)))
    ks=np.max(a,b)
    
    return ks, s, r


def lamb_normalize(cs_list):
    '''
    CS normalization following Lamb et al., 2006
    Input:
        cs_list: list, connectivity scores for drugs of interest
    Output:
        list, normalized connectivity scores
    '''
    
    p=max(cs_list)
    q=min(cs_list)
    
    def extr(x):
        if x>0:
            return p
        return -q
    
    return [x/extr(x) for x in cs_list]


def calculate_CS(a_up, a_down, b_up, b_down):
    '''
    Calculate connectivity score.
    from Lamb et al., 2006
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - CS: float, CS value in the interval [-1, 1]
    '''
    
    if a_up>b_up:
        ks_up=a_up
    else:
        ks_up=-b_up
        
    if a_down>b_down:
        ks_down=a_down
    else:
        ks_down=-b_down
    
    if ks_up*ks_down>0:
        return 0
    return ks_up-ks_down   


def calculate_CS_evil_twin(ks_up, ks_down):
    '''
    Calculate connectivity score. But its evil twin.
    For testing
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - CS: float, CS value in the interval [-1, 1]
    '''
    return ks_up-ks_down 


def calculate_RGES(a_up, a_down, b_up, b_down):
    '''
    Calculate Reverse Gene Expression Score (RGES).
    from Bin Chen et al., 2017
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - RGES: float, RGES value in the interval [-2, 2]
    '''
    
    if a_up>b_up:
        ks_up=a_up
    else:
        ks_up=-b_up
        
    if a_down>b_down:
        ks_down=a_down
    else:
        ks_down=-b_down
    
    return ks_up-ks_down   


def montecarlo_connectivity(s_up, s_down, r, n_iterations=1000, score_type='bin_chen'):
    '''
    Randomly sample RGES
    Input:
        - s_up: int, number of upregulated genes
        - s_down: int, number of downregulated genes
        - r: int, total number of genes in reference drug expression data
        - n_iterations: int, number of Monte Carlo iterations (default: 1000)
        - score_type: str, optional, which calculation to use. 
            options: ['bin_chen', 'sirota', 'lamb'], default: 'bin_chen'
    Output:
        - list of float, sampled RGES values
    '''
    
    random_RGES_list=[]
    drug_genes=np.arange(r)
    
    for i in range(n_iterations):
        
        # Generate random indexes for up and downregulated genes
        random_up_and_down_indexes=np.random.choice(r, s_up + s_down, replace=False)
         
        # Create random index dictionary:
        random_V_up=random_up_and_down_indexes[:s_up]
        random_V_down=random_up_and_down_indexes[-s_down:]
        
        # Compute random KS stats:
        random_a_up, random_b_up, _, _ = compute_KS(random_up_and_down_indexes[:s_up], random_V_up, drug_genes)
        random_a_down , random_b_down, _, _ = compute_KS(random_up_and_down_indexes[-s_down:], random_V_down, drug_genes)

        # calculate RGES evil twin:
        if score_type=='evil_twin':
            ks_up,  _, _ = compute_KS(random_up_and_down_indexes[:s_up], random_V_up, drug_genes)
            ks_down, _, _ = compute_KS(random_up_and_down_indexes[-s_down:], random_V_down, drug_genes)
            random_RGES_list.append(calculate_CS_evil_twin(ks_up, ks_down, random_b_up, random_b_down))
        
        # Calculate random RGES:
        if score_type=='bin_chen':
            random_RGES_list.append(calculate_RGES(random_a_up, random_a_down, random_b_up, random_b_down))
        if (score_type=='lamb') or(score_type=='sirota') :
            random_RGES_list.append(calculate_CS(random_a_up, random_a_down, random_b_up, random_b_down))

    if score_type=='lamb':
        random_RGES_list=lamb_normalize(random_RGES_list)

    return random_RGES_list

def random_disease_bin_chen_connectivity(drug_signature, rank_on='magnitude'):
    '''temp
    calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017, of a drug signature versus
   a random genes  ranking
    Input:
        - drug_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - rank_on: str ranking column. default='magnitude' 
                options: {'p_value', 'magnitude'}. Do not change unless you
                know what you're doing'
                    
    Output:
        -   measured_RGES: float, calculated RGES value
        -   p_value: float, p-value from Monte Carlo simulation
    '''
    
    
    if rank_on=='p_value':
        ranking_col_name=drug_signature.columns[2]
    if rank_on=='magnitude':
        ranking_col_name=drug_signature.columns[1]
    
    # generate random up and down regulated disease genes
    r= drug_signature.shape[0]
    l=np.random.randint(1,r)
    s_up = r-l
    s_down = l
    gene_id_list = drug_signature.gene_id
    shuffled_gene_ids = gene_id_list[np.random.choice(gene_id_list.index, len(gene_id_list), replace=False)]
    random_disease_up = shuffled_gene_ids[:s_up]
    random_disease_down = shuffled_gene_ids[-s_down:]
     
    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, pd.DataFrame(random_disease_up)  , ranking_col_name, 'gene_id')
    V_down = rank_genes(drug_signature, pd.DataFrame(random_disease_down), ranking_col_name, 'gene_id')
    
    # Compute random KS stats:
    a_up, b_up, _, _ = compute_KS(random_disease_up, V_up, drug_signature['gene_id'])
    a_down , b_down, _, _ = compute_KS(random_disease_down, V_down, drug_signature['gene_id'])
    
    # Compute RGES:
    measured_RGES_for_random_disease = calculate_RGES(a_up, a_down, b_up, b_down)
    
    return measured_RGES_for_random_disease

def bin_chen_connectivity(disease_signature, drug_signature, rank_on='magnitude'):
    '''calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017
    Input:
        - disease_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - drug_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - rank_on: str ranking column. default='magnitude' 
                options: {'p_value', 'magnitude'}. Do not change unless you
                know what you're doing'
                    
    Output:
        -   measured_RGES: float, calculated RGES value
        -   p_value: float, p-value from Monte Carlo simulation
    '''
    
    disease_p_val_col_name=disease_signature.columns[2]
    disease_signature_col_name=disease_signature.columns[1]
    
    if rank_on=='p_value':
        ranking_col_name=drug_signature.columns[2]
    if rank_on=='magnitude':
        ranking_col_name=drug_signature.columns[1]
    
    # # Get lists of up (down) regulated genes: older with 2vs
    disease_signature_up = disease_signature[disease_signature[disease_signature_col_name]>0]
    disease_signature_down = disease_signature[disease_signature[disease_signature_col_name]<0]
    
    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, disease_signature_up  , ranking_col_name, disease_p_val_col_name)  # CHECK THIS LINE
    V_down = rank_genes(drug_signature, disease_signature_down, ranking_col_name, disease_p_val_col_name)
    
    # Compute KS statistic for up and down regulated genes:
    a_up, b_up, s_up, r = compute_KS(disease_signature_up['gene_id'], V_up, drug_signature['gene_id'])
    a_down , b_down, s_down, r = compute_KS(disease_signature_down['gene_id'], V_down, drug_signature['gene_id'])
    
    # Compute RGES:
    measured_RGES = calculate_RGES(a_up, a_down, b_up, b_down)
    
    # Calculate two tailed p-value (no assumption on the direction of RGES 
    # between disease and drug) for measured RGES, using random sampling:
    n_iterations=1000 
    random_RGES_list = montecarlo_connectivity(s_up, s_down, r, n_iterations, score_type='bin_chen')
    p_value=np.sum(np.abs(np.array(random_RGES_list))>np.abs(measured_RGES))/n_iterations
    return measured_RGES, p_value

#%%     CALCULATE CONNECTIVITY FOR ONE DRUG DATA, FOR ONE PERTURBATION TIME
# see cs_batch for multiple drugs

if __name__=='__main__':
    
    def run_connectivity_score(DISEASE, mith, drug, pert_time ):
        '''
        wrapper function to load disease and drug data (for one drug), 
        and calculate their connectivity score.
        '''
        singature_uom = 'DE_log2_FC' if not mith else 'Perturbation'
        # load disease signature:
        disease_signature=load_disease_signature(DISEASE, mith=mith)
        disease_signature=disease_signature[['gene_id', singature_uom, 'adj.p.value']]
        drug_signature=load_single_drug_signature(drug, mith=mith)
        
        disease_common_index, drug_common_index = get_common_genes(disease_signature, drug_signature)  
        
        # filter disease signature with common genes:
        disease_signature=disease_signature.iloc[disease_common_index].reset_index(drop=True)
        print('common disease genes', len(disease_signature))
        # load drug signature
        p_val_str= 'adj.p.value' if not pert_time=='6h_24h' else 'p.value'
        columns_of_interest=['gene_id', singature_uom+'_'+pert_time, p_val_str+'_'+pert_time]
        cscore, p_value = bin_chen_connectivity(disease_signature, drug_signature[columns_of_interest], rank_on='magnitude')
        print(cscore, p_value)
        return cscore, p_value
    
    mith=True
    pert_time='6h_24h' 
    drug='ibuprofen' #'cortisone' nice value against ipf, ibuprofen, nice against als_NYGC
    print(DISEASE)
    print(drug)
    cscore, p_value = run_connectivity_score(DISEASE, mith, drug, pert_time)    
