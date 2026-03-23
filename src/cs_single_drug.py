# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:28:36 2026

#     CALCULATE CONNECTIVITY FOR ONE DRUG DATA, FOR ONE PERTURBATION TIME
#     for given disease 
#     config: 18 11 2025 IPF vs Catalano LINCS main pipeline
# see cs_batch for multiple drugs

@author: los4
"""

from conf import DISEASE
from loader import load_disease_signature, load_single_drug_signature
from connectivity_score import get_common_genes, bin_chen_connectivity

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