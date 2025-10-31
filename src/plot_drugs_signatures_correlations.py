#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:59:42 2025

@ author: L-F-S

Inspects correlations between DEG signatures and MITHrIL signatures 
for LINCS drug signatures:
    
- Creates and save images of the top 5 most and least correlated DEG vs mith 
drug signatures (and pvalues)

- Saves list of correlations between DEG and mith signatures and pvalues
     for all drugs in tsv in output/drug_signatures_correlations
"""

from scipy import stats
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from conf import DISEASE, IMG_DIR, MITH_IN_DRUG, MITH_OUT_DRUG,\
 TSR_OUT_DRUG,TSR_OUT_DISEASE, TSR_OUT_CSCORE, OOUT_DIR
from preprocessing_utils import get_drugs_list
from loader import load_disease_signature, load_single_drug_signature, load_drug_signatures
import time
import  matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# functions

def compare_deg_and_mith_genes(DEG_signatures, mith_signatures):
    '''
    Print number of common genes between mithril and DEG datasets, and return genes in common
    '''
    missing_genes=set.difference(set(DEG_signatures.gene_id),set(mith_signatures['gene_id']))
    print('missing genes from mithril', len(missing_genes))
    addional_mith_items=set.difference(set(mith_signatures['gene_id']),set(DEG_signatures.gene_id))
    print('addional_mith_items', len(addional_mith_items))
    common_genes=set.intersection(set(mith_signatures['gene_id']),set(DEG_signatures.gene_id))
    print('genes in common', len(common_genes))
    print('riprova: common genes+ additional mithril items == len(mithril_gene)?',len(common_genes)+len(addional_mith_items)==mith_signatures.shape[0])
    print('riprova: common genes+ missing mithril items == len(DEG_genees)?',len(common_genes)+len(missing_genes)==DEG_signatures.shape[0])

    # Get subset of signatures with common genes
    mith_common_genes_signatures=mith_signatures[mith_signatures['gene_id'].isin(common_genes)]
    DEG_common_genes_signatures=DEG_signatures[DEG_signatures.gene_id.isin(common_genes)]
    
    return DEG_common_genes_signatures, mith_common_genes_signatures


def prop_agreeing_signs(DEG_common_set, mith_common_set, sign_col):
    '''
    Print proportion of equal signs between DEG and mithril genes (only genes in common)
    '''
    return mith_common_set.set_index('gene_id').sort_index()[sign_col][mith_common_set.set_index('gene_id').sort_index()[sign_col]==DEG_common_set.set_index('gene_id').sort_index()[sign_col]].shape[0]/mith_common_set.shape[0]


def print_correlations(DEG_signature, mith_signature, \
                       DEG_col='DE_log2_FC', mith_col='Perturbation',data_name=DISEASE, scatter=True):
    '''
    plots scatterplot of DEG vs Mithril signatures for given variables (default DEG Fc vs Mithril Perturbation),
    and prints corresponding pearson and spearman coefficients
    
    Input signatures are assumed to have the same unique gene_ids
    '''
    
    mith_common_perturbations=mith_signature[mith_col]
    DEG_common_FC=DEG_signature[DEG_col]
    DEG_v_mith_pearson=stats.pearsonr(DEG_common_FC, mith_common_perturbations)
    DEG_v_mith_spearman=stats.spearmanr(DEG_common_FC,mith_common_perturbations)
    
    if scatter:
        '''Print correlations and plot scatter plot'''
        
        plt.scatter(DEG_common_FC,mith_common_perturbations)
        plt.xlabel('DEG '+DEG_col)
        plt.ylabel('mith '+mith_col)
        plt.title('FC vs mith signature for '+str(DEG_signature.shape[0])+' genes for '+data_name +' data')
        print('Pearson cc:',DEG_v_mith_pearson,'\nSpearman cc:', DEG_v_mith_spearman)
        
    return DEG_v_mith_pearson, DEG_v_mith_spearman


def plot_correlation(correlations_dataframe, correlation):
    plt.figure()
    correlations_dataframe[correlation].sort_values(ascending=False).plot(kind='bar')
    plt.xticks(rotation=75)
    plt.ylim(-1,1)
    plt.axhline(0, color='darkgrey', linestyle='--')
    plt.title(correlation)
    plt.savefig(IMG_DIR+correlation+'.pdf')
    return

def plot_10_most_correlated(correlations_dataframe, correlation):
    plt.figure()
    most_and_least_correlated=pd.concat([correlations_dataframe[correlation].sort_values(ascending=False).head(5), correlations_dataframe[correlation].sort_values(ascending=False).tail(5)])
    most_and_least_correlated.plot(kind='bar')
    plt.xticks(rotation=75)
    # plt.ylim(-1,1)
    plt.axhline(0, color='darkgrey', linestyle='--')
    plt.title(correlation)
    plt.savefig(IMG_DIR+correlation+'.pdf')
    return

def save_correlation(correlations_dataframe, correlation):
    if not os.isdir(OOUT_DIR+'drug_signatures_correlations'+os.sep):
        os.mkdir(OOUT_DIR+'drug_signatures_correlations'+os.sep)
    correlations_dataframe[correlation].sort_values(ascending=False).to_csv(OOUT_DIR+'drug_signatures_correlations/'+correlation+'.tsv', sep='\t', index_label='drug')
    return
    
    
#%%


if __name__=='__main__':
    
    
    # list of all Pearson and Spearman correlations
    
    pearsons_magnitude_6h=[]
    spearmans_magnitude_6h=[]
    
    pearsons_pvalue_6h=[]
    spearmans_pvalue_6h=[]
    
    
    pearsons_magnitude_24h=[]
    spearmans_magnitude_24h=[]
    
    pearsons_pvalue_24h=[]
    spearmans_pvalue_24h=[]
    
    
    pearsons_magnitude_6h_24h=[]
    spearmans_magnitude_6h_24h=[]
    
    pearsons_pvalue_6h_24h=[]
    spearmans_pvalue_6h_24h=[]
    
    
    # slice common genes only, from first drug. Since the DEG drugs
    # have the same genes, they will all have the same common genes with mith
    
    ex_drug_DEG_signature=load_single_drug_signature('ibuprofen', mith=False)
    ex_drug_mith_signature=load_single_drug_signature('ibuprofen', mith=True)
    
    
    
    drug_DEG_common, drug_mith_common = compare_deg_and_mith_genes(ex_drug_DEG_signature, ex_drug_mith_signature)  
    
    drug_DEG_common.sort_values(by='gene_id', inplace=True)
    drug_mith_common.sort_values(by='gene_id', inplace=True)
    drug_DEG_common_index=drug_DEG_common.index
    drug_mith_common_index=drug_mith_common.index
    
    
    drugs_list=get_drugs_list()
    


    
        
    for drug in drugs_list:
        
        print(drug)
        drug_DEG_common=load_single_drug_signature(drug,mith=False).iloc[drug_DEG_common_index]
        
        drug_mith_common=load_single_drug_signature(drug,mith=True).iloc[drug_mith_common_index]

        
        # calculating correlation coefficients with 6h signatures
        drug_pearson_magnitude, drug_spearman_magnitude = print_correlations(drug_DEG_common, drug_mith_common,\
                                                         DEG_col='DE_log2_FC_6h', mith_col='Perturbation_6h',data_name=drug,\
                                                         scatter=False)
        pearsons_magnitude_6h.append(drug_pearson_magnitude[0])
        spearmans_magnitude_6h.append(drug_spearman_magnitude[0])

    
        drug_pearson_pval, drug_spearman_pval = print_correlations(drug_DEG_common, drug_mith_common,\
                                                         DEG_col='p.value_6h', mith_col='p.value_6h',data_name=drug,\
                                                         scatter=False)
        pearsons_pvalue_6h.append(drug_pearson_pval[0])
        spearmans_pvalue_6h.append(drug_spearman_pval[0])
    
    
        # calculating correlation coefficients with 24h signatures
        drug_pearson_magnitude, drug_spearman_magnitude = print_correlations(drug_DEG_common, drug_mith_common,\
                                                         DEG_col='DE_log2_FC_24h', mith_col='Perturbation_24h',data_name=drug,\
                                                         scatter=False)
        pearsons_magnitude_24h.append(drug_pearson_magnitude[0])
        spearmans_magnitude_24h.append(drug_spearman_magnitude[0])
    
        drug_pearson_pval, drug_spearman_pval = print_correlations(drug_DEG_common, drug_mith_common,\
                                                         DEG_col='p.value_24h', mith_col='p.value_24h',data_name=drug,\
                                                         scatter=False)
        pearsons_pvalue_24h.append(drug_pearson_pval[0])
        spearmans_pvalue_24h.append(drug_spearman_pval[0])
    
        
        # calculating correlation coefficients with 6h_24h signatures
        drug_pearson_magnitude, drug_spearman_magnitude = print_correlations(drug_DEG_common, drug_mith_common,\
                                                         DEG_col='DE_log2_FC_6h_24h', mith_col='Perturbation_6h_24h',data_name=drug,\
                                                         scatter=False)
        pearsons_magnitude_6h_24h.append(drug_pearson_magnitude[0])
        spearmans_magnitude_6h_24h.append(drug_spearman_magnitude[0])
    
        drug_pearson_pval, drug_spearman_pval = print_correlations(drug_DEG_common, drug_mith_common,\
                                                         DEG_col='p.value_6h_24h', mith_col='p.value_6h_24h',data_name=drug,\
                                                         scatter=False)
        pearsons_pvalue_6h_24h.append(drug_pearson_pval[0])
        spearmans_pvalue_6h_24h.append(drug_spearman_pval[0])


    correlations_dataframe=pd.DataFrame(index=drugs_list)
    correlations_dataframe['pearsons_magnitude_6h']=pearsons_magnitude_6h
    correlations_dataframe['pearsons_magnitude_24h']=pearsons_magnitude_24h
    correlations_dataframe['pearsons_magnitude_6h_24h']=pearsons_magnitude_6h_24h
    
    correlations_dataframe['spearmans_magnitude_6h']=spearmans_magnitude_6h
    correlations_dataframe['spearmans_magnitude_24h']=spearmans_magnitude_24h
    correlations_dataframe['spearmans_magnitude_6h_24h']=spearmans_magnitude_6h_24h
    
    correlations_dataframe['pearsons_pvalue_6h']=pearsons_pvalue_6h
    correlations_dataframe['pearsons_pvalue_24h']=pearsons_pvalue_24h
    correlations_dataframe['pearsons_pvalue_6h_24h']=pearsons_pvalue_6h_24h
    
    correlations_dataframe['spearmans_pvalue_6h']=spearmans_pvalue_6h
    correlations_dataframe['spearmans_pvalue_24h']=spearmans_pvalue_24h
    correlations_dataframe['spearmans_pvalue_6h_24h']=spearmans_pvalue_6h_24h
    
    
    
    
    
    
#%%
    
if __name__=='__main__':   
    
    
    for correlation in correlations_dataframe.columns:
        print('plotting 10 least and most correlated signatures for', correlation)
        plot_10_most_correlated(correlations_dataframe, correlation)
        print('saving all signatures for', correlation)

        save_correlation(correlations_dataframe, correlation)