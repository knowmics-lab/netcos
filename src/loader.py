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
    
def load_landmark_gene_ids(path):
    landmark_genes = pd.read_csv(path / 'lincs_landmark.csv') 
    return list(landmark_genes.gene_id.sort_values())

###############################################################################
#
# Connectivity Score and Validation data
#
###############################################################################
def add_pert_id_to_cs( lincs_metadata_path,cs_df, cs_id_col='LINCS_id', metadata_id_col='id', metadata_pert_col='pert_id'):
    """
    Adds a 'pert_id' column to a CS dataframe by mapping LINCS_id via LINCS metadata.
    
    Parameters
    ----------
    cs_df : pd.DataFrame
        Connectivity score dataframe containing LINCS_id column
    lincs_metadata_path : str
        Path to lincs_sig_info_new.csv
    cs_id_col : str
        Column in cs_df (default 'LINCS_id')
    metadata_id_col : str
        Column in metadata corresponding to LINCS_id (default 'id')
    metadata_pert_col : str
        Column in metadata for pert_id (default 'pert_id')
    """

    if not metadata_pert_col in cs_df.columns:
        # load only what we need
        meta = pd.read_csv(lincs_metadata_path, usecols=[metadata_id_col, metadata_pert_col], dtype='str')
    
        # build mapping dict
        id_to_pert = dict(zip(meta[metadata_id_col], meta[metadata_pert_col]))
    
        # map
        cs_df[metadata_pert_col] = cs_df[cs_id_col].map(id_to_pert)
    
        # optional sanity check
        n_missing = cs_df['pert_id'].isna().sum()
        if n_missing > 0:
            print(f"Warning: {n_missing} LINCS_id values could not be mapped to pert_id")

    return cs_df

    
def load_drug_rankings(path, filename=None, pert_time='all',mith='mith', lincs_metadata_path =None):
    
    if not filename:
        filename=mith+'_connectivity_score.tsv'
    df= pd.read_csv(path/filename, sep='\t', header=0, dtype='str')

    if lincs_metadata_path:
        df = add_pert_id_to_cs(lincs_metadata_path, df) 
    
    if not pert_time == 'all':
        return df[df.perturbation_time==pert_time]
        
    return df

def translate_cl(cell_line):
    if cell_line == 'HT29':
        return 'HT-29'
    return cell_line

def load_IC50(ic50_file, cancer_type, cell_line, IC50_ONLY=True):#, median_IC50=False):
    
    log_dict = {}
    df=pd.read_excel(ic50_file,\
                     sheet_name=cancer_type, header =0, usecols=['pert_iname', 'pert_id','standard_value', 'standard_units',\
                                 'standard_type', 'cell_line','activity', 'standard_value_median'])
    #filter for cell line
    cell_line = translate_cl(cell_line)
    if cell_line is not None:
        log_dict['IC50_rows'] = len(df)
        df = df[df['cell_line'].astype(str).str.upper()==cell_line.upper()]
        log_dict['IC50_rows_filtered_by_cell'] = len(df)
    
    # IC50 filtering
    if IC50_ONLY:
        df=df[df.standard_type=='IC50']
        log_dict['IC50_rows_filtered_by_IC50'] = len(df)

    
    return df.sort_values(by='standard_value_median'), log_dict