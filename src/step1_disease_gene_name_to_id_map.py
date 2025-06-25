#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 11:13:41 2025

@ author: L-F-S

Convert Disease fold change gene names into
gene IDs in the tsr output,
and creates MITHRiL input
"""
import os

import pandas as pd
import pickle
import numpy as np
import collections
from conf import DISEASE,MITH_IN_DISEASE, TSR_OUT_DISEASE, DICT_DIR, alias_2geneid

######################
# DISEASE DATA
#####################
# preprocessing to remove duplicates

disease_gene_signature_data=pd.read_csv(TSR_OUT_DISEASE+os.sep+DISEASE+'_signature.csv', sep=';', decimal=',',header=0)
#%% remove duplicates
genes_not_in_alias=[]
disease_genes_id_to_symbol_univocous={}
duplicated_ids_with_different_symbols=collections.defaultdict(list)
for symbol in disease_gene_signature_data['gene']:
    try:
        gene_id=str(alias_2geneid[symbol])
#        if not gene_id in disease_genes_id_to_symbol_univocous.keys():
        disease_genes_id_to_symbol_univocous[gene_id] = symbol #actually not univocous yet, there are duplicates
        duplicated_ids_with_different_symbols[gene_id].append(symbol)
    except:

        genes_not_in_alias.append(symbol)
print(len(genes_not_in_alias), 'gene symbols without gene id ') 

#only keep duplicates:
topop=[]
for gene_id, symbols  in duplicated_ids_with_different_symbols.items():
    if len(symbols)==1:
        topop.append(gene_id)
for gene_id in topop:
    duplicated_ids_with_different_symbols.pop(gene_id)
print(len(duplicated_ids_with_different_symbols.keys()), 'duplicated ids')

#%%

def print_duplicated(genes,df):
    FCs=[]
    for gene_symbol in genes:
        fc=disease_gene_signature_data[disease_gene_signature_data['gene']==gene_symbol]['DE_log2_FC'].iloc[0]
#        np.round(fc,3)
#        print(fc, type(fc))
        FCs.append(np.round(fc,3))
    if not len(np.unique(np.array(FCs)))==1:
#        print('duplicated:',genes)

        return 1
    return 0

count_duplicated_genes_with_different_FC=0
for gene_id, duped_symbols in duplicated_ids_with_different_symbols.items():
    count_duplicated_genes_with_different_FC+=print_duplicated(duped_symbols, disease_gene_signature_data)
print('tot duplicated genes with idfferent FC', count_duplicated_genes_with_different_FC)
genes_id_to_symbol_univocous={gene_id:disease_genes_id_to_symbol_univocous[gene_id] for gene_id in topop}

#%%save dictionary for later translation into gene names if needed:
filename=DICT_DIR+DISEASE+'_gene_id_to_symbol.pkl'
with open(filename, 'wb') as f:
    pickle.dump(genes_id_to_symbol_univocous, f)
    
# filename=DICT_DIR+DISEASE+'_gene_id_to_symbol.pkl'
# with open(filename, 'rb') as f:
#     genes_id_to_symbol_univocous=pickle.load(f)

#%% Create gene id mithril input, by keeping only relevant genes
gene_id_mith_input=[]
for gene_id in topop: # slow
    fc=disease_gene_signature_data[disease_gene_signature_data['gene']==disease_genes_id_to_symbol_univocous[gene_id]]['DE_log2_FC'].iloc[0]
    gene_id_mith_input.append(gene_id+'\t'+str(fc))

#%% save mithril input
    
f=open(MITH_IN_DISEASE+DISEASE+'_signature_gene_id.mi','w')
f.write(('\n').join(gene_id_mith_input))
f.close()

#%% Convert the original gene signature data

# remove gene symbols not in mapping
disease_gene_signature_data = disease_gene_signature_data[~disease_gene_signature_data.gene.isin(genes_not_in_alias)]

# Remove duplicates
disease_gene_signature_data['gene_id']=disease_gene_signature_data['gene'].apply(lambda gene : str(alias_2geneid[gene]))
disease_gene_signature_data = disease_gene_signature_data[disease_gene_signature_data.gene_id.isin(topop)]

# Save file
disease_gene_signature_data.to_csv(TSR_OUT_DISEASE+DISEASE+os.sep+DISEASE+'_signature_gene_id.csv', sep=';',decimal=',', index=None)

