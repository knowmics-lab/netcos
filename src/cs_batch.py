# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 12:08:22 2025

@author: L-F-S

Parallel calculation of conectivity score and other similarity measures
for a list of drugs vs a disease
"""

import os
import sys
import numpy as np
import pandas as pd
import time
from scipy import stats
from connectivity_score import get_common_genes, calc_connectivity_score_with
from loader import load_disease_signature, load_single_signature_cs_input,\
    load_landmark_gene_ids
from preprocessing_utils import get_signature_ids_list_from_cs_input
from conf import DISEASE, CS_OUT, CS_IN_DRUG, CS_IN_DISEASE, cs_batch_threads,\
    disease_run_name, cs_on_LM, LINCS_BC_DATA, cs_mith, LOGS_DIR,\
        cell_line_run_name, cs_log_filename, lincs_metadata_path,\
            CS_ON_PATHWAYS,CS_METHOD, make_cs_filename
from preprocessing_utils import get_chunk_indexes
from logger import append_run_metadata
from joblib import Parallel, delayed

from pathlib import Path
from datetime import datetime
import socket

# def append_cs_run_log(metadata_path,row_dict):
#     """
#     Append one run row to cs_runs.tsv, creating the file if it does not exist.
#     """
#     metadata_path = Path(metadata_path)
#     metadata_path.parent.mkdir(parents=True, exist_ok=True)

#     row_df = pd.DataFrame([row_dict])

#     if metadata_path.exists():
#         old_df = pd.read_csv(metadata_path, sep='\t')
#         out_df = pd.concat([old_df, row_df], ignore_index=True)
#     else:
#         out_df = row_df

#     out_df.to_csv(metadata_path, sep='\t', index=False)

###############################################################################
#             CALCULATE CONNECTIVITY FOR DEG/MITH DATA
###############################################################################

def add_pert_id_to_cs( lincs_metadata_path,cs_df,
                      cs_id_col='LINCS_id',
                      metadata_id_col='id',
                      metadata_pert_col='pert_id'):
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
        n_missing = cs_df[metadata_pert_col].isna().sum()
        if n_missing > 0:
            print(f"Warning: {n_missing} LINCS_id values could not be mapped to pert_id")

    return cs_df

def run_connectivity_score_drugs_batch(disease_run_name, mith, drugs_list, n, i1, i2, \
                                       rank_on='magnitude', save_file=False, cs_on_LM=False,\
                                        cs_on_pathways=False):
    '''
    wrapper function to load disease and drug data (for a  batch of drugs in drugs_list[i1:i2]), 
    and calculate their connectivity score, for given drugs and LINCS drug perturbation times.
    Also calculates Pearson correlation
    coefficient, Spearman correlation coefficient, and cosine similarity, between
    disease and drug.
    Output: a dataframe of connectivity scores and  pvalues, for all combinations of conditions
    '''
    start=time.time()
    run_stats = {}
    
    if not cs_on_pathways:
        singature_uom = 'DE_log2_FC' if not mith else 'Perturbation'
        id_col = 'gene_id'
        pval_col = 'adj.p.value'
        columns_of_interest = [id_col, singature_uom, pval_col]
    else:
        singature_uom = "Corrected Accumulator"
        id_col = 'ID'
        pval_col = 'Adjusted pValue'
        columns_of_interest = [id_col, singature_uom, pval_col]

    # load disease signature:
    disease_signature=load_disease_signature(disease_run_name, mith=mith, pathway=CS_ON_PATHWAYS)
    
    #log variable
    run_stats['n_disease_before_common'] = len(disease_signature)
    disease_signature=disease_signature[columns_of_interest]
    
    run_stats["disease_value_col"] = singature_uom
    run_stats["disease_pval_col"] = pval_col
    run_stats["disease_id_col"] = id_col
    run_stats["drug_value_col"] = singature_uom
    run_stats["drug_pval_col"] = pval_col
    run_stats["drug_id_col"] = id_col

    # Discard disease genes that are 
    # not found in drug genes
    # DO NOT simply select common genes, as this would impact
    # connectrivity calculations
    # load one drug signature (all durg signatures have the same genes in the same order)
    
    drug_signature= load_single_signature_cs_input(drugs_list[0], CS_IN_DRUG, pathway=CS_ON_PATHWAYS)
    # log variable
    run_stats['n_drug_before_common'] = len(drug_signature)
    
    # log variables
    run_stats['n_disease_after_second_lm'] = None
    run_stats['n_drug_after_second_lm'] = None
    
    if cs_on_LM:
        lm_gene_ids = [str(x) for x in load_landmark_gene_ids(LINCS_BC_DATA)]
        disease_signature = disease_signature[disease_signature[id_col].isin(lm_gene_ids)].reset_index(drop=True)
        drug_signature = drug_signature[drug_signature[id_col].isin(lm_gene_ids)].reset_index(drop=True)
        run_stats['n_disease_after_second_lm'] = len(disease_signature)
        run_stats['n_drug_after_second_lm'] = len(drug_signature)
    
    # filter common genes
    disease_common_index, drug_common_index = get_common_genes(disease_signature, drug_signature, id_col=id_col)     
    # filter disease signature with common genes:
    # Note: indexing every time with pandas .iloc is slower than writing
    # and loading binaries for filtered dataframes, for each drug
    # however this is not relevant, as we only need to filter
    # disease data, while we can take all drug data.
    disease_signature=disease_signature.iloc[disease_common_index].reset_index(drop=True)
    
    #log variables
    run_stats['n_disease_after_common'] = len(disease_signature)
    run_stats['n_drug_after_common'] = len(drug_common_index)

    
    # Initialize dataframe:
    data = []
    
    
    # Calculate connectivity score between disease and drugs:     
    for drug in drugs_list[i1:i2]:
        
        # load drug signature
        
        drug_signature=load_single_signature_cs_input(drug, CS_IN_DRUG, pathway=CS_ON_PATHWAYS)
        if cs_on_LM:
            drug_signature = drug_signature[drug_signature['gene_id'].isin(lm_gene_ids)].reset_index(drop=True)
        
            
        # Calculate connectivity score between drug and disease:
        cscore, p_value = calc_connectivity_score_with[CS_METHOD](
            disease_signature, \
            drug_signature[columns_of_interest], rank_on=rank_on, id_col=id_col)
        
        # Calculate other correlations (between common genes)
        disease_signature_value=disease_signature[singature_uom]
        drug_signature_value=drug_signature[singature_uom].loc[drug_common_index]
        
        pearson=stats.pearsonr(disease_signature_value, drug_signature_value)
        spearman=stats.spearmanr(disease_signature_value, drug_signature_value)
        cosine_sim=np.dot(np.array(disease_signature_value),np.array(drug_signature_value))/(np.linalg.norm(np.array(disease_signature_value))*np.linalg.norm(np.array(drug_signature_value)))
        
        # carry metadata through if present
        row = {
            "disease": disease_run_name,
            "LINCS_id":drug,
            "connectivity_score": cscore,
            "cs_p_value": p_value,
            "pearson": pearson[0],
            "pearson_p_value": pearson[1],
            "spearman": spearman[0],
            "spearman_p_value": spearman[1],
            "cos_sim": cosine_sim,
        }
        for col in ["pert_iname", "pert_desc", "pert_id", "sig_id", "cell_id", "pert_time", "pert_dose", "pert_type", "is_gold", "DrugBank.ID"]:
            if col in drug_signature.columns:
                row[col] = drug_signature[col].iloc[0]
        
        data.append(row)
    
        # data.append([DISEASE, drug, cscore, p_value, pearson[0], pearson[1], spearman[0], spearman[1], cosine_sim]) # other correlations, genes subset
    # connectivity_data= pd.DataFrame(data, columns=["disease","pert_id","connectivity_score","cs_p_value",'pearson','pearson_p_value','spearman','spearman_p_value','cos_sim']) #add other correlations, genes subset
    connectivity_data= pd.DataFrame(data)
    elapsed = time.time() - start
    print('total elapsed time for batch',str(n),'/',n_jobs,' of', len(drugs_list[i1:i2]), 'drugs: ', elapsed)
    
    run_stats['elapsed']=elapsed
    run_stats['n_drugs_in_batch']=len(drugs_list[i1:i2])
    
    # save connectivity scores
    if save_file:
        
        if not os.path.exists(CS_OUT):
            os.mkdir(CS_OUT)
        
        if not cs_on_pathways:
            connectivity_dataset_filename= str(i1)+'_'+str(i2)+'_DEG_connectivity_score.tsv' if not mith else  CS_OUT+str(i1)+'_'+str(i2)+'_mith_connectivity_score.tsv'
        
        else:
            connectivity_dataset_filename= CS_OUT+str(i1)+'_'+str(i2)+'_mith_pw_connectivity_score.tsv'
        connectivity_data.to_csv(CS_OUT/connectivity_dataset_filename, sep='\t', index=False)
    
    print('total elapsed time for batch of', len(drugs_list[i1:i2]),' drugs: ', time.time()-start)
    return connectivity_data, run_stats
#%% Parallel run for LINCS data
import itertools
if __name__=="__main__":
    
    #Hyperparameter master
    
    cs_on_LMs = [0, 1]
    # cs_miths = [0, 1]
    # CS_ON_PATHWAYSs = [False]
    CS_METHODs = ['bin_chen', 'bin_chen_disease_sorted']
    
    
    #############################
    # BEGIN CS CALCULATION FOR SINGLE SET OF PARAMETERS FOR ALL DRUGS
    for cs_on_LM, CS_METHOD in itertools.product(cs_on_LMs, CS_METHODs):
        print("cs on LM:", str(cs_on_LM))
        print("cs method:", str(CS_METHOD))
        connectivity_dataset_filename, cs_id, = make_cs_filename(cs_mith, cs_on_LM, CS_ON_PATHWAYS, CS_METHOD)
        
        start_total = time.time()
        # Set calculation for MITHrIL data
        mith=cs_mith
        drugs_list=get_signature_ids_list_from_cs_input(CS_IN_DRUG)
        lm_flag = cs_on_LM
        
        # set n of cores for parallel execution
        n_jobs= cs_batch_threads # int(sys.argv[1]) #16
        if n_jobs>len(drugs_list):
            raise ValueError('n_jobs:',n_jobs,' > len(drugs_list):',len(drugs_list),'! Reduce size of n_jobs')
        
        # get chunk size
        chunk_size=int(len(drugs_list)/n_jobs)
        last_chunk_size=len(drugs_list)%n_jobs
        parallel_indexes=[]
        for i in range(n_jobs):
            i1, i2 =  get_chunk_indexes(i, chunk_size)
            parallel_indexes.append((i1,i2))
        
        print('n jobs', n_jobs)
        print('len drug list', len(drugs_list))
        print('chunk size',chunk_size)
        print('last chunk size', last_chunk_size)
        print('mith tag:', mith, '\ndisease:', disease_run_name)
        print('pathway based signature:', CS_ON_PATHWAYS)
        
        # results = run_connectivity_score_drugs_batch(disease_run_name, mith, drugs_list, i1, i2, cs_on_LM=lm_flag)
        # cs_df=results[0]
        # first_stats=results[1]
        results = Parallel(n_jobs=n_jobs)(delayed(run_connectivity_score_drugs_batch)\
                                     (disease_run_name, mith, drugs_list, n, i1, i2, cs_on_LM=lm_flag, cs_on_pathways=CS_ON_PATHWAYS)\
                                for n, (i1,i2) in enumerate(parallel_indexes))
        
        if last_chunk_size>0:
            last_batch=run_connectivity_score_drugs_batch(disease_run_name, mith, drugs_list, (n_jobs+1), i2,len(drugs_list), cs_on_LM=lm_flag, cs_on_pathways=CS_ON_PATHWAYS)
            results.append(last_batch)
        
        # build dataframe    
        cs_df=pd.concat([x[0] for x in results], ignore_index=True)
        
        # take stats from first batch as representative of filtering sizes
        first_stats = results[0][1]
        
        #To test:
        cs_df = add_pert_id_to_cs(lincs_metadata_path, cs_df) 
        
        #
        total_elapsed = time.time() - start_total
        
        # write results
        if not os.path.exists(CS_OUT):
                    os.mkdir(CS_OUT)
        #%% write file
        
        # now =  datetime.now()
        # datetime_string = now.strftime("%d_%m_%Y_%H_%M")
        # cs_filename = datetime_string+'_DEG_connectivity_score.tsv' if not cs_mith else  datetime_string+'_mith_connectivity_score.tsv'
        # connectivity_dataset_filename=CS_OUT/cs_filename
    
        
        cs_df.to_csv(connectivity_dataset_filename, sep='\t', index=False)
        
        #%% write run log
        
        
        metadata_row = {
            "cs_run_id": cs_id,
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "hostname": socket.gethostname(),
            "disease_run_id": disease_run_name,
            "drug_run_id": cell_line_run_name,
            "mith": int(mith),
            "cs_on_LM": int(lm_flag),
            "CS_ON_PATHWAYS": CS_ON_PATHWAYS,
            "rank_on": "magnitude",\
            "drug_id_col": first_stats["drug_id_col"] ,
            "disease_id_col": first_stats["disease_id_col"] ,
            "drug_value_col": first_stats["drug_value_col"],
            "drug_pval_col": first_stats["drug_pval_col"],
            "disease_value_col": first_stats["disease_value_col"],
            "disease_pval_col": first_stats["disease_pval_col"] ,
            "CS_METHOD": CS_METHOD,
            "n_drugs_total": len(drugs_list),
            "n_results_rows": len(cs_df),
            "n_disease_genes_before_common": first_stats["n_disease_before_common"],
            "n_drug_genes_before_common": first_stats["n_drug_before_common"],
            "n_disease_genes_after_second_lm": first_stats["n_disease_after_second_lm"],
            "n_drug_genes_after_second_lm": first_stats["n_drug_after_second_lm"],
            "n_disease_genes_after_common": first_stats["n_disease_after_common"],
            "n_drug_genes_after_common": first_stats["n_drug_after_common"],
            "n_jobs": n_jobs,
            "elapsed_sec_total": total_elapsed,
            "drug_input_dir": str(CS_IN_DRUG),
            "disease_input_dir": str(CS_IN_DISEASE),
            "output_file": str(connectivity_dataset_filename),
        }
    
        append_run_metadata(cs_log_filename, metadata_row)
            
        #############################
        # END CS CALCULATION FOR SINGLE SET OF PARAMETERS FOR ALL DRUGS