#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:41:12 2024

@ author: L-F-S
"""

#%%
from conf import MITH_OUT_DRUG
import os
import numpy as np
import subprocess
from preprocessing_utils import get_drugs_list



#%% Get all drugs
drugs_list=get_drugs_list()

# set n of cores
cores=30

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
    subprocess.Popen(['python','preprocessing_utils.py',str(i1), str(i2)])

subprocess.Popen(['python','preprocessing_utils.py',str(len(drugs_list)-chunk_size),str(len(drugs_list))])

