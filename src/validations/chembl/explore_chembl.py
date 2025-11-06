# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:33:55 2025

@author: L-F-S
"""
#%%

from conf import cell_lines_chembl, DATA_DIR
from validations.chembl.chembl_loader import load_raw_chembl
import pandas as pd
import os
import numpy as np

cell_line = cell_lines_chembl[1]
act_df= load_raw_chembl(cell_line, activity_index=True)


#%% get InChIKey 
# map molecule_chembl_id to InChIKey
formulas = act_df.canonical_smiles.unique()
ex_molecule = formulas[0]


#%% Load LINCS data for given drugs for given cell line
import h5py
LINCS_file = DATA_DIR+'LINCS-GSE92742'+os.sep+'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
#%%
def explore_hdf5(name, obj):
    print(name)

with h5py.File(LINCS_file, "r") as f:
    f.visititems(explore_hdf5)
    
    # 0
    # 0/DATA
    # 0/DATA/0
    # 0/DATA/0/matrix
    # 0/META
    # 0/META/COL
    # 0/META/COL/id
    # 0/META/ROW
    # 0/META/ROW/id
#%%
with h5py.File(LINCS_file, "r") as f:
    # Explore the internal structure
    print(list(f.keys()))  # Typically includes '0' or 'DATA'
    print(f.name)
    
    data = f['/0/DATA/0/matrix']
    
    print(list(f['0/META/COL'].keys()))
    # Read sample IDs (row IDs in the matrix)
    # A sample is a single expression profile 
    # for one cell line, exposed to one drug (or control)
    # at a given dose, for one time point.
    sample_ids = f['0/META/COL/id'][:].astype(str)
    # data = f['0']
    print(data.shape)  # e.g., (1319138, 12328)
    chunk = data[:1000, :1000]  # Get first 1000 columns (samples)
#%%

inst_info_file = DATA_DIR + 'LINCS-GSE92742' + os.sep + 'GSE92742_Broad_LINCS_inst_info.txt'

inst_df = pd.read_csv(inst_info_file, sep='\t', low_memory=False)
#%%
def filter_index_for_cell_Line(inst_df, cell_line):
    # Filter for your target cell lines
    filtered_inst_df = inst_df[inst_df['cell_id'].isin(cell_lines)]
    print(f"Found {len(filtered_inst_df)} matching samples.")
    target_ids = filtered_inst_df['inst_id'].values
    # Create mapping: inst_id → index in matrix
    id_to_index = {id_: i for i, id_ in enumerate(sample_ids)}
    
    # Get row indices for selected samples
    target_indices = [id_to_index[inst_id] for inst_id in target_ids if inst_id in id_to_index]
    print(f"Matched {len(target_indices)} indices out of {len(target_ids)} requested.")
    
    # sort indices (for h5py)
    target_indices = np.array(target_indices)
    sorted_idx = np.argsort(target_indices)
    sorted_indices = target_indices[sorted_idx]
    return target_indices, sorted_indices, sorted_idx, filtered_inst_df

cell_line = cell_lines[0]
target_indices, sorted_indices,sorted_idx, filtered_inst_df = filter_index_for_cell_Line(inst_df, cell_line)
#%%  Only extract data for given filters
with h5py.File(LINCS_file, "r") as f:
    matrix = f['0/DATA/0/matrix']
    selected_data = matrix[sorted_indices, :] 
    print(selected_data.shape)

unsorted_data = selected_data[np.argsort(sorted_idx)]
#%%
with h5py.File(LINCS_file, "r") as f:
    gene_ids = f['0/META/ROW/id'][:].astype(str)
expr_df = pd.DataFrame(unsorted_data, index=target_indices, columns=gene_ids)

#%%
filtered_inst_df = filtered_inst_df.set_index('inst_id')
# Group sample IDs by pert_id and pert_dose
grouped = filtered_inst_df.groupby(['pert_id', 'pert_dose']).apply(lambda df: df.index.tolist())
#%% Example: print groups and number of samples
for (drug, dose), samples in grouped.items():
    print(f"Drug: {drug}, Dose: {dose}, Samples: {len(samples)}")

#%%