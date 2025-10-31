# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:33:55 2025

@author: L-F-S
"""

from conf import CHEMBL_INPUT_DATA_DIR, cell_lines
import pandas as pd

def load_raw_chembl(cell_line, activity_index = False):
    index_col = 'activity_id' if activity_index else False
    return pd.read_csv(CHEMBL_INPUT_DATA_DIR+cell_line+'_activity.tsv', sep='\t', index_col=index_col)

#%%
