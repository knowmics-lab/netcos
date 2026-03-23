#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:24:39 2025

@ author: L-F-S
"""
import os

import pandas as pd
import pickle
from conf import CS_DIR, TSR_OUT_DRUG, TSR_OUT_DISEASE, CS_IN_DISEASE, CS_IN_DRUG
from preprocessing_utils import get_drugs_list

def load_disease_signature(DISEASE, mith=False):
    '''
    inputs:
        DISEASE: str disease symbol (directory name)

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.

    returns a pandas.Dataframe where rows are genes.
    '''
    
    if not mith:
        filename=os.path.join(TSR_OUT_DISEASE,DISEASE+'_signature_gene_id.csv')
        return pd.read_csv(filename, sep=';',decimal=',', dtype={'gene_id':'str'})
   
    else:
        filename=os.path.join(CS_IN_DISEASE,DISEASE+'_mith3_signature.csv')
        return pd.read_csv(filename, sep='\t')

def load_single_drug_signature(drug, mith=False, pkl=True):
    '''
    Loads DEG drug wise drug signature (with unique gene id rows)

    inputs:
        DISEASE: str disease symbol (directory name)

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.

    returns a pandas.Dataframe where rows are genes.
    '''
    if not mith:
        
        filename=os.path.join(TSR_OUT_DRUG,'LINCS','metanalysis_drug_wise_filtered',drug+'_metanalysis')
        if not pkl:
            return pd.read_csv(filename+'.csv', sep='\t',  dtype={'gene_id':'str'})
        else:
            with open(filename+'.pkl', 'rb') as f:
                data=pickle.load(f)
            return data
                
    else:
        
        filename=os.path.join(CS_IN_DRUG,drug+'_metanalysis')
        
        if not pkl:
            return pd.read_csv(filename+'.csv', sep='\t')
        else:
            with open(filename+'.pkl', 'rb') as f:
                data=pickle.load(f)
            return data
           

def load_drug_signatures(mith=False, pkl=True):
    '''
    inputs:

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.
    returns a pandas.Dataframe which is a concatenation by row
    of all drug signature data. rows are genes signature for given drug 
    (genes will appear multiple time, once for every drug)
    '''
    
    print('loading all drug signatures...')
    
    drugs_list=get_drugs_list(mith)
    
    if not mith:
        DEG_drug_list=[]
        for n, drug_file in enumerate(drugs_list): 
            drug=drug_file.split('_')[0]
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    else:
        DEG_drug_list=[]
        for n, drug_file in enumerate(drugs_list):
            drug=drug_file.split('_')[0]
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    print(n+1, 'drugs loaded')
    return pd.concat(DEG_drug_list)

def load_single_signature_cs_input(signature_id, cs_input_dir):
    filename = os.path.join(cs_input_dir, f"{signature_id}.pkl")
    with open(filename, "rb") as f:
        return pickle.load(f)