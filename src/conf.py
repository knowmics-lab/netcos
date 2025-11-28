#!'+os.sep+'usr'+os.sep+'bin'+os.sep+'env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

18-11-2025 conf file for chembl valdiation pipeline.
for HEPG 24h
"""

import os
import pickle
from local import BASE_DIR
####################
# FLAGS
####################
 
# Disease symbol. This will be used in the pipeline to 
# identify the disease. It will be in file names and directories
# !! IMPORTANT: keep spaces around '=' and spaces before cmments for pipeline to work  
DISEASE = 'temp'#'LIHC' #'ipf'  # 'als_NYGC' # 
cell_line = 'MCF7'
pert_time = '24h'


################
# DIRECTORIES
################
import local 
# only override if it exists in local.py
if hasattr(local, "DATA_DIR"):
    DATA_DIR = local.DATA_DIR
else:
    DATA_DIR=BASE_DIR+'data'+os.sep
    
LINCS_DIR=BASE_DIR+'data'+os.sep+'LINCS-GSE92742'+os.sep
DICT_DIR=DATA_DIR+'dictionaries'+os.sep

MITH_APP="/home/signorini/mithril3/app-3.0.0-SNAPSHOT.jar"

MITH_DIR=BASE_DIR+'MITHrIL'+os.sep
MITH_OUT=MITH_DIR+'output'+os.sep
MITH_OUT_DISEASE=MITH_DIR+'output'+os.sep+'disease_signature_2025'+os.sep
MITH_OUT_DRUG=MITH_DIR+'output'+os.sep+'drug_signature_2025'+os.sep+cell_line+'_'+pert_time+os.sep

MITH_IN=MITH_DIR+'input'+os.sep
MITH_IN_DISEASE=MITH_DIR+'input'+os.sep+'disease_signature'+os.sep
MITH_IN_DRUG=MITH_DIR+'input'+os.sep+'drug_signature'+os.sep


TSR_DIR=BASE_DIR+'tsr'+os.sep
TSR_OUT=TSR_DIR+'output'+os.sep
TSR_OUT_DISEASE=TSR_OUT+'disease_signature'+os.sep+DISEASE+os.sep
TSR_OUT_DRUG=TSR_OUT+'drug_signature'+os.sep
TSR_OUT_DRUG_META=TSR_OUT_DRUG+'LINCS_lorenzo'+os.sep+'metanalysis_mith3_drug_wise'+os.sep
TSR_OUT_CSCORE=TSR_OUT+'connectivity_score'+os.sep

# Connectivity score dirs
CS_DIR=BASE_DIR+'connectivity_score'+os.sep
CS_IN_DRUG=CS_DIR+'input'+os.sep+'drug_signature_2025'+os.sep+cell_line+'_'+pert_time+os.sep
CS_IN_DISEASE=CS_DIR+'input'+os.sep+'disease_signature_2025'+os.sep+DISEASE+os.sep
CS_OUT=CS_DIR+'output'+os.sep+DISEASE+'_2025'+os.sep

# other outputs
IMG_DIR=BASE_DIR+'imgs'+os.sep
OOUT_DIR=BASE_DIR+'other_outputs'+os.sep


##########################################################################
# Drug signature calculations using LINCS1000 dataset
##########################################################################


# LINCS1000 data filename (warning: big file)
LINCS_file = LINCS_DIR+os.sep+'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
#LINCS1000 metadata filename:
inst_info_file = LINCS_DIR+os.sep+ 'GSE92742_Broad_LINCS_inst_info.txt'
#LINCS1000 gene info filename:
GENE_INFO_FILE = LINCS_DIR+os.sep+ 'GSE92742_Broad_LINCS_gene_info.txt'
# landmark genes filename:
BING_GENES = LINCS_DIR+os.sep+"bing_gene_symbols.csv"  # optional

########################################################################
# MITHrIL parameters
########################################################################
mith_threads = 10
mith_organism = 'hsa_2025'


###############################################################################
# Chembl validation parameters
###############################################################################
VAL_DIR = BASE_DIR+os.sep+'validations'+os.sep 
# Chembl directory
CHEMBL_BASE_DIR = VAL_DIR+'chembl'+os.sep
CHEMBL_INPUT_DATA_DIR = CHEMBL_BASE_DIR+'chembl_input'+os.sep
cell_lines_chembl = ["MCF7", "HepG2", "HT29"]
#cell_lines_to_tumors = {"MCF7":, "HepG2":'LIHC', "HT29"}
# other validations
# SaveRunner directory
SR_DIR = VAL_DIR + os.sep + 'saverunner_comparison' + os.sep

####################
# Utility functions
####################

# load map of gene symols : gene ids
# more symbols may map to the same gene id

# using signorini's mapping (see network_signing algorithm)
alias_2geneid_filename=DICT_DIR+'alias_2geneid.pkl'
with open(alias_2geneid_filename, 'rb') as f:
    alias_2geneid = pickle.load(f)

def map_name_to_id(genename):
    if genename in alias_2geneid.keys():
        return alias_2geneid[genename]
    return genename

# # using alaimo's mapping
# symbol_to_id_filename=BASE_DIR+'other_data'+os.sep+'symbol_2geneid.pkl'
# with open(symbol_to_id_filename, 'rb') as f:
#     symbol_2geneid = pickle.load(f)

# def map_name_to_id(genename):
#     if genename in symbol_2geneid.keys():
#         return symbol_2geneid[genename]
#     return genename
