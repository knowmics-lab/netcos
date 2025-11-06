# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 17:12:15 2025

@author: L-F-S
"""

from conf import cell_lines_chembl,  alias_2geneid
import pandas as pd
import os
import numpy as np

def map_assay_id_and_matrix_idx(inst_df, sample_ids):
    id_to_index = {id_: i for i, id_ in enumerate(sample_ids)}
    inst_df.set_index('inst_id', inplace = True)
    inst_df["assay_index"] = inst_df.index.map(id_to_index)
    return id_to_index, inst_df

def filter_by_min_cell_lines_per_drug(inst_df, min_cells=5):
    """
    Keep only rows corresponding to drugs (pert_id) that have been tested
    in at least `min_cells` distinct cell lines.
    
    Parameters:
        inst_df (pd.DataFrame): Metadata dataframe with 'pert_id' and 'cell_id' columns.
        min_cells (int): Minimum number of distinct cell lines per drug.

    Returns:
        pd.DataFrame: Filtered dataframe.
    """

    # Count distinct cell lines per drug
    cell_counts = inst_df.groupby('pert_id')['cell_id'].nunique()
    
    # Keep only drugs tested in at least `min_cells` cell lines
    valid_drugs = cell_counts[cell_counts >= min_cells].index
    filtered_df = inst_df[inst_df['pert_id'].isin(valid_drugs)].copy()

    print(f"Remaining rows: {len(filtered_df)}")
    return filtered_df



def filter_by_cell_line(inst_df, cell_lines=[]):
    if len(cell_lines)==0:
        return inst_df
    return inst_df[inst_df['cell_id'].isin(cell_lines)]

def filter_by_pert_time(inst_df, allowed_times=[], only_int = True):
    """
    Keep only pert_id + pert_dose groups that have samples for *all* allowed_times.
    """
    
    if len(allowed_times)==0:
        return inst_df
    
    # Ensure we're only considering allowed times first
    df = inst_df[inst_df['pert_time'].isin(allowed_times)].copy()
    
    if only_int:
        # AGGIUNGERE QUESTA COSA COSI COSATA??
        # AGGIUNGERE LA COSA CHE UNO HA
        # CHE USANDO POCHE CELL LINES, SE NON FILTRO PER INTERSEIONE , HO 3.7K
        # DRUGS, SE FILTRO , NE HO TIPO 122, PERCHE
        # PROIBABILMENTE LUI FILTRAVA COSI PENSANDO A ANALISI SU TUTTE LE CELL LINES
        # MENTRE QUI NO:
        # Count unique times per (pert_id, pert_dose)
        print('taking intersection of perturbation times')
        grouped = df.groupby(['pert_id', 'pert_dose'])['pert_time'].nunique()
        # Keep only groups with all required time points
        valid_groups = grouped[grouped == len(allowed_times)].index
        # Filter inst_df by valid groups
        mask = df.set_index(['pert_id', 'pert_dose']).index.isin(valid_groups)
        df = df[mask]
    return df

def filter_by_dose(inst_df, doses=[]):
    if len(doses)==0:
        return inst_df
    return inst_df[inst_df['pert_dose'].isin(doses)]

def main_filter(inst_df, id_to_index, min_cell_lines=5, cell_lines=cell_lines_chembl,
                pert_times=[6, 24], doses=[10], only_int=True):

    # add/remove filters here:

    print(f"Filtering by {len(cell_lines)} cell lines...")
    inst_df = filter_by_cell_line(inst_df, cell_lines)
    
    print(f"Filtering by perturbation times: {pert_times}")
    inst_df = filter_by_pert_time(inst_df, pert_times, only_int)
    
    print(f"Filtering by doses: {doses}")
    inst_df = filter_by_dose(inst_df, doses)
    
    print(f"Filtered drugs to those used in ≥{min_cell_lines} cell lines.")
    inst_df = filter_by_min_cell_lines_per_drug(inst_df, min_cells=min_cell_lines)
    
    # TEMPPPP
    # print('TEMPORARILY slicing only first rows, for computability')
    # inst_df = inst_df.iloc[:100]
    #########
    
    print(f"Found {len(inst_df)} matching samples.")
    # Drop missing assay index (shouldn’t happen if map is complete)
    inst_df = inst_df.dropna(subset=['assay_index'])
    inst_df['assay_index'] = inst_df['assay_index'].astype(int)
    print(f"Found {len(inst_df)} matching samples.")
    print("Found", len(inst_df['pert_id'].unique()), 'perturbagen (drug) ids')

    
    # Sort for h5py slicing (must be ascending)
    sorted_idx = np.argsort(inst_df['assay_index'].values)
    sorted_indices = inst_df['assay_index'].values[sorted_idx]

    return inst_df, sorted_indices, sorted_idx

# TODO add genes preprocessing functions from LINCS preprocessing in chembl..