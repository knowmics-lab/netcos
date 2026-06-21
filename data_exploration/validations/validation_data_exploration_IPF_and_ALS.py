# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 11:08:32 2025

@author: L-F-S
"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('../modules')
import pandas as pd
from conf import DISEASE, BASE_DIR, CS_OUT
#%% search for drugs 
# IPF DRUGS:
    # Nintedanib 
    # pirfenidone but hey we dont know how effective they are do we?

CS_OUT = BASE_DIR + 'connectivity_score\\output'+os.sep+'ipf_2025'+os.sep # change disease here directly
mith_cs_data=pd.read_csv(CS_OUT+'mith_connectivity_score.tsv', sep='\t', index_col='drug')
mith_cs_data.loc['pirfenidone', ['perturbation_time','connectivity_score', 'cs_p_value' ]]
mith_cs_data.loc['nintedanib', ['perturbation_time','connectivity_score', 'cs_p_value' ]]

# ALS_NYGC DRUGS:
    # riluzole
    # edaravone/89-25-8/ Radicut/ 3-METHYL-1-PHENYL-2-PYRAZOLIN-5-ONE/Norphenazone   NOT found in dataset
CS_OUT = BASE_DIR + 'connectivity_score\\output'+os.sep+'als_NYGC_2025'+os.sep # change disease here directly
mith_cs_data=pd.read_csv(CS_OUT+'mith_connectivity_score.tsv', sep='\t', index_col='drug')
mith_cs_data.loc['riluzole', ['perturbation_time','connectivity_score', 'cs_p_value' ]]
# mith_cs_data.loc['edaravone', ['perturbation_time','connectivity_score', 'cs_p_value' ]]
