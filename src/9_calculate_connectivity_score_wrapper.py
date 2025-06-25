#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:41:12 2024

@ author: L-F-S
"""

#%%
from conf import CS_OUT
import pandas as pd
import subprocess
from preprocessing_utils import get_drugs_list
import os


#%% Get all drugs
drugs_list=get_drugs_list()

# set mithril flag (0:DEG data, 1:mithril data)
mith=0

# set n of cores
cores=16

# set ranking metric
rank_on='magnitude'

#%%

# get chunk size
chunk_size=int(len(drugs_list)/cores)
last_chunk_size=len(drugs_list)%cores

def get_chunk_indexes(i,chunk_size):
    i1=chunk_size*i
    i2=chunk_size*i+chunk_size
    return i1, i2

#%% 
for i in range(cores):
    i1,i2=get_chunk_indexes(i, chunk_size)
    subprocess.Popen(['python','connectivity_score.py',str(i1), str(i2),str(mith),rank_on])

subprocess.Popen(['python','connectivity_score.py',str(len(drugs_list)-chunk_size),str(len(drugs_list)), str(mith), rank_on])

