# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 10:52:55 2025

@author: los4
"""


import os

import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from preprocessing_utils import get_drugs_list
from conf import DISEASE, MITH_IN_DRUG, MITH_OUT_DRUG, CS_OUT, \
 TSR_OUT_DRUG,TSR_OUT_DISEASE, TSR_OUT_CSCORE,IMG_DIR
from loader import load_disease_signature, load_single_drug_signature, load_drug_signatures
from plot_drugs_signatures_correlations import compare_deg_and_mith_genes,  prop_agreeing_signs, print_correlations, plot_correlation, plot_10_most_correlated
print(DISEASE)

mith_cs_data_mg=pd.read_csv(CS_OUT+'mith_connectivity_score_magnitude.tsv', sep='\t')
DEG_cs_data_mg=pd.read_csv(CS_OUT+'DEG_connectivity_score_magnitude.tsv', sep='\t')
mith_cs_data=pd.read_csv(CS_OUT+'mith_connectivity_score.tsv', sep='\t')
DEG_cs_data=pd.read_csv(CS_OUT+'DEG_connectivity_score.tsv', sep='\t')

#%%
DEG_cs_data_mg.drop_duplicates(subset=['drug','perturbation_time'], inplace=True)
DEG_cs_data.drop_duplicates(subset=['drug','perturbation_time'], inplace=True)
mith_cs_data_mg.drop_duplicates(subset=['drug','perturbation_time'], inplace=True)
mith_cs_data.drop_duplicates(subset=['drug','perturbation_time'], inplace=True)

#%%
DEG_cs_data_mg.to_csv(CS_OUT+'DEG_connectivity_score_magnitude.tsv', sep='\t', index=False)
DEG_cs_data.to_csv(CS_OUT+'DEG_connectivity_score.tsv', sep='\t', index=False)
mith_cs_data_mg.to_csv(CS_OUT+'mith_connectivity_score_magnitude.tsv', sep='\t', index=False)
mith_cs_data.to_csv(CS_OUT+'mith_connectivity_score.tsv', sep='\t', index=False)
