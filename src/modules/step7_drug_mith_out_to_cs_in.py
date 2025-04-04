#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

Maps  Mithril 3 batch output into
tsr metanalysis-like input matrix.
"""

import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import pickle
import numpy as np
import pandas as pd
import time
from conf import MITH_OUT_DRUG, CS_IN_DRUG, TSR_OUT_DRUG
from preprocessing_utils import get_drugs_list


def remove_special_characters(df):
    special_chars=';/:*?\"|'
    for char in special_chars:
        indexes=df[df['Gene Name'].str.contains('\\'+char)].index
        if len(indexes)>0:
            df.drop(index=indexes, inplace=True)
    return df


#%%

def mith_out_to_cs_in(drugs_list,i1=None,i2=None, save_csv=False):
    print('Mapping mith3 output into tsr connectivity input: \n\
              Removing pathway duplicates, Merging three timepoint datasets in a single file and filtering ')
    start=time.time()
    tot_drugs=len(drugs_list)
    for h, drug in enumerate(drugs_list[i1:i2]):
        
        output_filename_csv=TSR_OUT_DRUG+'LINCS/metanalysis_mith3_drug_wise/'+drug+'_metanalysis.csv'   
        output_filename_pkl=CS_IN_DRUG+drug+'_metanalysis.pkl'   
        if not os.path.isfile(output_filename_pkl): 
        
        
            drug_start=time.time()
            print(drug, h+1, 'of', tot_drugs)
            #6h
            mith_perturb_signature_file_6h=drug+'_6h.perturbations.txt'
            mith_perturb_signature_6h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h, sep='\t', index_col=False, engine='python', usecols=['Gene Id', 'Gene Name', 'Perturbation', 'pValue', 'adj_pValue'])
            mith_perturb_signature_6h.drop_duplicates(inplace=True, keep='first') # Duplicate gene ids for pathways (reduces row number from 250k to 14k)
            mith_perturb_signature_6h.rename(columns={'Gene Id':'gene_id', 'Perturbation':'Perturbation_6h', 'pValue':'p.value_6h', 'adj_pValue':"adj.p.value_6h"}, inplace=True)
            mith_perturb_signature_6h['drug'] = drug
            mith_perturb_signature_6h['t.value_like_statistic_6h']=mith_perturb_signature_6h['p.value_6h'].apply(lambda p_value : 2001*(1-p_value))
        
            #24h
            mith_perturb_signature_file_24h=drug+'_24h.perturbations.txt'
            mith_perturb_signature_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_24h, sep='\t', index_col=False, engine='python', usecols=['Gene Id', 'Perturbation', 'pValue', 'adj_pValue'])
            mith_perturb_signature_24h.drop_duplicates(inplace=True, keep='first')
            mith_perturb_signature_24h.rename(columns={'Gene Id':'gene_id', 'Perturbation':'Perturbation_24h', 'pValue':'p.value_24h', 'adj_pValue':"adj.p.value_24h"}, inplace=True)
            mith_perturb_signature_24h['t.value_like_statistic_24h']=mith_perturb_signature_24h['p.value_24h'].apply(lambda p_value : 2001*(1-p_value))
        
            mith_perturb_signature_meta=mith_perturb_signature_6h.merge(mith_perturb_signature_24h, on='gene_id', how='left')
           
            #6h 24h
            mith_perturb_signature_file_6h_24h=drug+'_6h_24h.perturbations.txt'
            mith_perturb_signature_6h_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h_24h, sep='\t', index_col=False, engine='python', usecols=['Gene Id','Perturbation', 'pValue', 'adj_pValue'])
            mith_perturb_signature_6h_24h.drop_duplicates(inplace=True, keep='first')
            mith_perturb_signature_6h_24h.rename(columns={'Gene Id':'gene_id', 'Perturbation':'Perturbation_6h_24h', 'pValue':'p.value_6h_24h', 'adj_pValue':"adj.p.value_6h_24h"}, inplace=True)
            mith_perturb_signature_6h_24h['t.value_like_statistic_6h_24h']=mith_perturb_signature_6h_24h['p.value_6h_24h'].apply(lambda p_value : 2001*(1-p_value))
            
            mith_perturb_signature_meta=mith_perturb_signature_meta.merge(mith_perturb_signature_6h_24h, on='gene_id', how='left')
        
            mith_perturb_signature_meta=remove_special_characters(mith_perturb_signature_meta)
            print('drug', drug, 'processed in', np.round(time.time()-drug_start, 1))
            
            print('writing mit3 connectivity input metanalysis files for drug', drug)
            if save_csv:
                # to save a human-readable .csv file. 
                # will slow down the computation and use more space
                mith_perturb_signature_meta.to_csv(output_filename_csv, sep='\t', index=False)
            with open(output_filename_pkl, 'wb') as f:
                pickle.dump(mith_perturb_signature_meta, f)
    print('time elapsed:', time.time()-start)
    return


#%%

if __name__=='__main__':
    from joblib import Parallel, delayed
    
    # Get all drugs
    drugs_list=get_drugs_list()
    
    # set n of cores
    cores=2
    
    # get chunk size
    chunk_size=int(len(drugs_list)/cores)
    last_chunk_size=len(drugs_list)%cores
    
    def get_chunk_indexes(i,chunk_size):
        i1=chunk_size*i
        i2=chunk_size*i+chunk_size
        return i1, i2
    
    parallel_indexes=[]
    for i in range(cores):
        i1,i2=get_chunk_indexes(i, chunk_size)
        parallel_indexes.append((i1,i2))
     
    
    
    # i1=int(sys.argv[1])
    # i2=int(sys.argv[2])
    save_csv=False#int(sys.argv[3]) # change to True if you want to save csv readable files
    # todo fallo carino che puoi chiamarlo wrapper se lo vuoi fare
    # parallelo, e invece se no lo fai iterativo
    # e elimina il wrapper 6_mapblablawrapper, xke lo farai sul notebook
    # e vedi se eliminare il step 8 cioe vedi se il pkl creato da questo
    # mapmith3 aggiornato a cui ha iaggiunto il pickl, funziona lanciando 
    # il connectivity score. bel lavoro pronto per il pomeriggio
    # mith_out_to_cs_in(drugs_list,i1=0,i2=10, save_csv=False)
    results = Parallel(n_jobs=2)(delayed(mith_out_to_cs_in)\
                                 (drugs_list, i1, i2,save_csv=False)\
                            for i1,i2 in parallel_indexes)