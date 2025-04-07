#!'+os.sep+'usr'+os.sep+'bin'+os.sep+'env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S
"""

import os
import pickle

####################
# FLAGS
####################
 
# Disease symbol. This will be used in the pipeline to 
# identify the disease. It will be in file names and directories
# !! IMPORTANT: keep spaces around '=' and spaces before cmments for pipeline to work  
DISEASE = 'ipf'  # 'als_NYGC' # 


################
# DIRECTORIES
################

# BASE_DIR=os.sep+'home'+os.sep+'signorini'+os.sep+'drug_repurposing'+os.sep
# BASE_DIR=os.sep+'home'+os.sep+'lorenzo'+os.sep+'unict_2024-25'+os.sep+'drug_repurposing'+os.sep
BASE_DIR='G:'+os.sep+'Il mio Drive'+os.sep+'unict 2024-25'+os.sep+'drug_repurposing'+os.sep+''
# BASE_DIR=os.sep+'home'+os.sep+'NetCos'+os.sep+'modules'+os.sep

DATA_DIR=BASE_DIR+'data'+os.sep
DICT_DIR=DATA_DIR+'dictionaries'+os.sep

MITH_APP="/home/signorini/mithril3/app-3.0.0-SNAPSHOT.jar"

MITH_DIR=BASE_DIR+'MITHrIL'+os.sep
MITH_OUT=MITH_DIR+'output'+os.sep
MITH_OUT_DISEASE=MITH_DIR+'output'+os.sep+'disease_signature'+os.sep
MITH_OUT_DRUG=MITH_DIR+'output'+os.sep+'drug_signature'+os.sep

MITH_IN=MITH_DIR+'input'+os.sep
MITH_IN_DISEASE=MITH_DIR+'input'+os.sep+'disease_signature'+os.sep
MITH_IN_DRUG=MITH_DIR+'input'+os.sep+'drug_signature'+os.sep


TSR_DIR=BASE_DIR+'tsr'+os.sep
TSR_OUT=TSR_DIR+'output'+os.sep
TSR_OUT_DISEASE=TSR_OUT+'disease_signature'+os.sep+DISEASE+os.sep
TSR_OUT_DRUG=TSR_OUT+'drug_signature'+os.sep
TSR_OUT_DRUG_META=TSR_OUT_DRUG+'LINCS_lorenzo'+os.sep+'metanalysis_mith3_drug_wise'+os.sep
TSR_OUT_CSCORE=TSR_OUT+'connectivity_score'+os.sep

# Joined model dirs
IMG_DIR=BASE_DIR+'imgs'+os.sep
CS_DIR=BASE_DIR+'connectivity_score'+os.sep
CS_IN_DRUG=CS_DIR+'input'+os.sep+'drug_signature'+os.sep
CS_IN_DISEASE=CS_DIR+'input'+os.sep+'disease_signature'+os.sep+DISEASE+os.sep
CS_OUT=CS_DIR+'output'+os.sep+DISEASE+os.sep

# landmark genes file:
    
landmark_genes_filepath=TSR_DIR+'data'+os.sep+'LINCS-GSE92742'+os.sep+'landmark_genes.csv'



####################
# Utility functions
####################

# load map of gene symols : gene ids
# more symbols may map to the same gene id

# using signorini's mapping (see network_signing algorithm)
alias_2geneid_filename=DICT_DIR+'alias_2geneid.pkl'
with open(alias_2geneid_filename, 'rb') as f:
    alias_2geneid = pickle.load(f)

# # using alaimo's mapping
# symbol_to_id_filename=BASE_DIR+'other_data'+os.sep+'symbol_2geneid.pkl'
# with open(symbol_to_id_filename, 'rb') as f:
#     symbol_2geneid = pickle.load(f)


# def map_name_to_id(genename):
#     if genename in symbol_2geneid.keys():
#         return symbol_2geneid[genename]
#     return genename


