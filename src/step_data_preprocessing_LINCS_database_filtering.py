# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 16:17:13 2025

@author: L-F-S
"""

from conf import cell_lines, DATA_DIR, alias_2geneid
from tsr_python.LINCS_preprocessing import map_assay_id_and_matrix_idx, main_filter
import pandas as pd
import os
import numpy as np
import h5py

# Paramedters (to move to conf.py file)
LINCS_file = DATA_DIR+'LINCS-GSE92742'+os.sep+'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
inst_info_file = DATA_DIR + 'LINCS-GSE92742' + os.sep + 'GSE92742_Broad_LINCS_inst_info.txt'


#%% load sample ids 
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
    
## FILTERING TESTS:
#%% filtering
# Catalano's filtering, for 3 cell lines
from conf import cell_lines
cell_lines = [cell_line.upper() for cell_line in cell_lines]
pert_times=[6, 24]
doses=[10]
min_cell_lines=0 

filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 100249 matching samples.
# Found 5787 pert_ids
len(filtered_inst_df.pert_iname.unique())
# 5301 pert_inames
#%% no filtering (should return all data). It does

cell_lines = []
pert_times=[]
doses=[]
min_cell_lines=0
 
filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 1319138 matching samples.
# Found 43526 pert_ids
#%% 
cell_lines = ["HEPG2"]
pert_times=[6,24]
doses=[10]
min_cell_lines=0
 
filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 1393 matching samples.
# Found 122 pert_ids
#%%

cell_lines = []
pert_times=[6,24]
doses=[10]
min_cell_lines=5 

filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 370258 matching samples.
# Found 5155 pert_ids
len(filtered_inst_df.pert_iname.unique())
# but 4844 pert_inames..

# and if I change min_cell_lines filtering to work on pert_iname instead of pert_id 
# I get a bit more even, which makes sense because there
# are less pert_inames...
#%%  
metadata_cols = ['rna_plate', 'rna_well', 'pert_id', 'pert_iname', 'pert_type',
       'pert_dose', 'pert_dose_unit', 'pert_time', 'pert_time_unit', 'cell_id',
       'assay_index']
  
pert_types=['ctl_vehicle', 'trt_cp', 'trt_lig', 'ctl_vector', 'trt_sh',
       'ctl_untrt', 'trt_oe', 'trt_oe.mut']

inst_df['pert_id'][inst_df['pert_type']=='ctl_vehicle'].unique()
# array(['DMSO', 'PBS', 'H2O'])
inst_df['pert_id'][inst_df['pert_type']=='ctl_vector'].unique()
# quindi c'e anche PBS...
# ma