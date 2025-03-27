"""
Created on Wed Nov 13 09:55:46 2024

@ author: L-F-S

extract unique gene-wise perturbation values from MITHrIL output,
for disease data
"""

import pandas as pd
from conf import  MITH_OUT_DISEASE, CS_IN_DISEASE, DISEASE   


def mith_to_tsr(mith_data):
    '''map mithril output to tsr output'''
    
    mith_data.drop_duplicates('Gene Id', keep='first', inplace=True)
    
    tsr_like=mith_data[['Gene Id','Gene Name','Perturbation','pValue', 'adj_pValue']].rename(columns={'pValue':'p.value','Gene Name':'gene','adj_pValue' :'adj.p.value', 'Gene Id':'gene_id'})

    
    return tsr_like


#%% Disease mithril3
# Open MITHRIL condition OutPut
mith_file=MITH_OUT_DISEASE+DISEASE+'_mith3.perturbations.txt'
mith_data=pd.read_csv(mith_file, sep='\t')


#%%

tsr_like=mith_to_tsr(mith_data)

#write file
tsr_output_name=DISEASE+'_mith3_signature.csv'
mith_condition_filename=TSR_OUT_DISEASE+DISEASE+'/'+tsr_output_name
tsr_like.to_csv(mith_condition_filename,sep=';', header=True, index=False)


