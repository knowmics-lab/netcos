# -*- coding: utf-8 -*-
"""
Created on Mon May 12 12:53:20 2025

@author: los4
"""


import os

import pandas as pd
from conf import DISEASE, BASE_DIR, CS_OUT

def filter_top_drugs(correlations_dataframe, pert_time='6h_24h',\
                                   correlation='connectivity_score',\
                                   p_value_col='cs_p_value', p_value=0.01,\
                                   connectivity_score=-1.95    ):
    
    '''
    correlations_dataframe is a dataframe with as index the label of interest (=drugs)
    correlation: str, correlation metric to use. Default: connectivity_score.
                                                other options: (pearson,
                                                                spearman,cos_sim).
    pert_time: str defines the experimental measure time point from initial
                    perturbation of medium. Default: 6h_24h (metanalysis data)
                                            other options: (6h, 24h)
    p_value=pvalue threshold to consider, if pvalue exists. default: 0.01)
    '''
    correlations_dataframe=correlations_dataframe[correlations_dataframe.perturbation_time==pert_time]
    correlations_dataframe=correlations_dataframe[[correlation,p_value_col]]
    return correlations_dataframe[correlations_dataframe[correlation]<connectivity_score][correlations_dataframe[p_value_col]<p_value].sort_values(by=correlation)

#%%
disease='ipf_2025'# als_NYGC' #+'_2025'
if __name__=='__main__':
    CS_OUT = BASE_DIR + 'connectivity_score\\output\\'\
        +os.sep+disease+os.sep # change disease here directly
    print(CS_OUT)
    print(DISEASE)
    mith_cs_data=pd.read_csv(CS_OUT+'mith_connectivity_score.tsv', sep='\t', index_col='drug')
    top_drugs=filter_top_drugs(mith_cs_data, p_value=0.05)
    print(top_drugs)

#%% find a single drug
drug='empagliflozin'

print(mith_cs_data.loc[mith_cs_data.index == drug,['perturbation_time', 'connectivity_score', 'cs_p_value']])
