#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:38:29 2025

@ author: L-F-S

create a version of a file list of bing_genes for the input of a common gene list between drug and disease
by using the genes of them ithril output which, for some reaso are not the same as the inèut 
(see log for details, it s not just metabolites and rnas added, but also genes removed)

"""
import os

import sys
import numpy as np
import pandas as pd
import pickle
from conf import BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG, TSR_OUT_DRUG, TSR_DIR

f=open(BASE_DIR+'other_data/mith_out_genes_list.txt')
genes=[x.strip() for x in f.readlines()]
f.close()
#%%
data=pd.DataFrame(genes, columns=['pr_gene_symbol'])
#%%
f=open(TSR_DIR+'/data/LINCS-GSE92742/landmark_genes.csv')
lmg=[x.strip() for x in f.readlines()]
f.close()
lmg=lmg[1:]
def is_landmark(g):
    if g in lmg:
        return 1
    return 0
data['pr_is_lm']=data['pr_gene_symbol'].apply(lambda g : is_landmark(g))
#%%
f=open(TSR_DIR+'/data/LINCS-GSE92742/bing_gene_symbols.csv')
bg=[x.strip() for x in f.readlines()]
f.close()
bg=bg[1:]
def is_bing(g):
    if g in bg:
        return 1
    return 0
data['pr_is_bing']=data['pr_gene_symbol'].apply(lambda g : is_bing(g))
#%%

data.reset_index(inplace=True)
data.rename(columns={'index':'pr_gene_id'}, inplace=True)
data['pr_gene_id']=data.index #placeholdes
data

data.to_csv(TSR_DIR+'/data/LINCS-GSE92742/mith_out_genes_list.txt', sep=';', header=True,index=False)
