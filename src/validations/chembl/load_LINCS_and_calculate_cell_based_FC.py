# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:33:55 2025

@author: L-F-S
"""

from conf import cell_lines_chembl, DATA_DIR, alias_2geneid
import pandas as pd
import os
import numpy as np
import h5py
from tsr_python.LINCS_preprocessing import map_assay_id_and_matrix_idx, main_filter

# Paramedters (to move to conf.py file)
LINCS_file = DATA_DIR+'LINCS-GSE92742'+os.sep+'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
inst_info_file = DATA_DIR + 'LINCS-GSE92742' + os.sep + 'GSE92742_Broad_LINCS_inst_info.txt'

cell_lines = [cell_line.upper() for cell_line in cell_lines_chembl]
print('temp using only one cell line for computability')
cell_lines = [cell_lines[1]] # 'HEPG2' has the least amount of data, so easier to debug temporary
pert_times=[6, 24]
doses=[10]
#%%

# def explore_hdf5(name, obj):
#     print(name)

# with h5py.File(LINCS_file, "r") as f:
#     f.visititems(explore_hdf5)
    
#     # 0
#     # 0/DATA
#     # 0/DATA/0
#     # 0/DATA/0/matrix
#     # 0/META
#     # 0/META/COL
#     # 0/META/COL/id
#     # 0/META/ROW
#     # 0/META/ROW/id


#%% temp
filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, cell_lines=cell_lines , pert_times=pert_times, doses=doses, only_int=True)

#%% Load LINCS data for given drugs for given cell line

with h5py.File(LINCS_file, "r") as f:
    # Explore the internal structure
    print(list(f.keys()))  # Typically includes '0' or 'DATA'
    print(f.name)
    
    matrix = f['/0/DATA/0/matrix']
    
    print(list(f['0/META/COL'].keys()))
    # Read sample IDs (row IDs in the matrix)
    # A sample is a single expression profile 
    # for one cell line, exposed to one drug (or control)
    # at a given dose, for one time point.
    sample_ids = f['0/META/COL/id'][:].astype(str)
    
    gene_ids = f['0/META/ROW/id'][:].astype(str)
    print(matrix.shape)  # (1319138, 12328)
        
    # Create mapping between assay id and index in matrix
    print('Loading metadata')
    inst_df = pd.read_csv(inst_info_file, sep='\t', low_memory=False)
    id_to_index, inst_df = map_assay_id_and_matrix_idx(inst_df, sample_ids)
    
    # Filter assay matrix by rows
    filtered_inst_df, sorted_indices, sorted_idx = \
        main_filter(inst_df, id_to_index,  cell_lines=cell_lines, pert_times=pert_times, doses=doses)

    #  Only extract data for given filters (h5py datasets only
    # allow for one dimension filtering at a time)
    # Warning: might be a big file
    print('extracting row indices from matrx...')
    selected_rows = matrix[sorted_indices,:] 
    print('Matrix of shape', selected_rows.shape)

#%% Filter by given genes
print('Filtering by bing genes')
bing_ids = pd.read_csv(DATA_DIR + 'LINCS-GSE92742' + os.sep + 'bing_genes.csv', header=0, sep=';', usecols=[0], dtype=np.str_).pr_gene_id.to_list()
genes_in_bing_index = np.where(np.isin(gene_ids,bing_ids))[0]
selected_data = selected_rows[:,genes_in_bing_index]

# ATTENZIONE: questo step e' diverso da Catalano!! commentato fuori
# quindi non lo metto
# landmark_genes = pd.read_csv(DATA_DIR + 'LINCS-GSE92742' + os.sep + 'landmark_genes.csv', header=0, sep=';', usecols=[0], dtype=np.str_).gene.to_list()
# landmark_ids = [alias_2geneid[gene] for gene in landmark_genes]
# filtered_gene_ids = set(genes_in_bing_index).union(set(landmark_ids))
# genes_in_bing_adn_landmark_index = np.where(np.isin(bing_ids,landmark_ids))[0]
# selected_data = selected_rows[:,genes_in_bing_adn_landmark_index]

del selected_rows

#restore data to 
# unsorted_data = selected_data[np.argsort(sorted_idx)]

expr_df = pd.DataFrame(selected_data, index=sorted_indices, columns=genes_in_bing_index)
print('filtered for bing genes, matrix of shape', expr_df.shape )

    
#%% get control assays
control_candidates = [pid for pid in filtered_inst_df['pert_id'].unique() if 'control' in pid.lower() or 'dmso' in pid.lower()]
print(control_candidates)
control_mask = filtered_inst_df['pert_id'] == 'DMSO'
control_samples = filtered_inst_df[control_mask]

print(f"Number of control samples: {len(control_samples)}")