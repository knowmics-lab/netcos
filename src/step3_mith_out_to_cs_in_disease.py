"""
Created on Wed Nov 13 09:55:46 2024

@ author: L-F-S

extract unique gene-wise perturbation values from MITHrIL output,
for disease data, and convert into connectivity score calculation input
"""

import pandas as pd
import os


def mith_to_cs_in(mith_file):
    '''map mithril output to tsr output'''
    
    # Read MITHrIL perturbation output, taking into account MITHrIL bug in 
    # header: there are 8 columns, but 7 column names. last column is
    # adjusted p value, but the column name is missing
    mith_data=pd.read_csv(mith_file, sep='\t', header=None, skiprows=1)
    mith_data.drop_duplicates(2, keep='first', inplace=True)
    
    cs_in=mith_data[[2,3,4,6,7]].rename(columns={2:'gene_id',3:'gene',4:'Perturbation',6:'p.value',7:'adj.p.value'})
    
    return cs_in


def write_cs_input(cs_data, cs_input_name, CS_IN_DISEASE):
    if not os.path.exists(CS_IN_DISEASE):
        os.mkdir(CS_IN_DISEASE)
    cs_data.to_csv(CS_IN_DISEASE+cs_input_name,sep='\t', header=True, index=False)
    print('Connectivity score disease input file for MITHrIL perturbation data written at',
          CS_IN_DISEASE,'\nfile name:',cs_input_name)
    return

#%%
if __name__=='__main__':
    
    from conf import  MITH_OUT_DISEASE, CS_IN_DISEASE, DISEASE   
    
    print(DISEASE)
    # read MITHrIL perturbation output
    mith_file=MITH_OUT_DISEASE+DISEASE+'_mith3.perturbations.txt'
    cs_data=mith_to_cs_in(mith_file)
    
    # write CS input for MITHrIL data
    cs_input_name=DISEASE+'_mith3_signature.csv'
    write_cs_input(cs_data, cs_input_name, CS_IN_DISEASE)
    

