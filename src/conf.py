#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

18-11-2025 conf file for chembl valdiation pipeline.
for HEPG cell line, LIHC disease, using landmark genes from BinChen2017 supll data
"""

from pathlib import Path
import pickle
from local import BASE_DIR
from datetime import datetime
# !! IMPORTANT: keep spaces around '=' and spaces before cmments for pipeline to work  


##############################################################################
# LINCS1000 level 4 data parameters
##############################################################################

landmark_disease = False  # only select landmark genes from signature data 
LM_flag_disease = ''
if landmark_disease:
    LM_flag_disease = '_LM'

# is_gold = 1 # only select high-quality perturbagens. Default:1 for all versions (33k perturbagens x 2 time steps = 66k perturbagens)
# pert_time = '6h'
landmark_drug = False  # only select landmark genes from signature data (irrelevant for Bin Chen data since they are already landmark only)
LM_flag_drug = ''
if landmark_drug:
    LM_flag_drug = '_LM'
cell_line = 'HT29'#'all'#
cell_line_run_name = cell_line+LM_flag_drug
mith_input_file = 'LINCS' +cell_line+LM_flag_drug+'.mi'



####################
# Disease data parameters
####################
 
# Disease symbol. This will be used in the pipeline to 
# identify the disease. It will be in file names and directories
# name of cell line for chembl validation
diseases_of = {'HEPG2':'LIHC',
               'MCF7':'BRCA',
               'HT29':'COAD',
               'all':'any'}

DISEASE = diseases_of[cell_line] #"ipf"#

disease_run_name = DISEASE + LM_flag_disease

###############################################################################
# DIRECTORIES and files
###############################################################################
import local 
# only override if it exists in local.py
if hasattr(local, "DATA_DIR"):
    DATA_DIR = local.DATA_DIR
else:
    DATA_DIR=BASE_DIR / 'data'
    
DICT_DIR=DATA_DIR / 'dictionaries'
LOGS_DIR = BASE_DIR / 'logs' 


# MITHrIL Directories

MITH_APP="/home/signorini/mithril3/app-3.0.0-SNAPSHOT.jar"

MITH_DIR=BASE_DIR / 'MITHrIL'
MITH_OUT=MITH_DIR / 'output'
MITH_OUT_DISEASE=MITH_DIR / 'output' / 'disease_signature_2025'
MITH_OUT_DRUG=MITH_DIR / 'output' / 'drug_signature_2025' / Path(cell_line_run_name) #Path(cell_line+'_'+pert_time)

MITH_IN=MITH_DIR / 'input'
MITH_IN_DISEASE=MITH_DIR / 'input' / 'disease_signature'
MITH_IN_DRUG=MITH_DIR / 'input' / 'drug_signature'

##########################################################################
# Bin Chen 2017 data directories
##########################################################################

if hasattr(local, "BC_DATA"):
    BC_DATA = local.BC_DATA
else:
    BC_DATA = DATA_DIR / 'BinChen2017'

if hasattr(local, "LINCS_BC_DATA"):
    LINCS_BC_DATA = local.LINCS_BC_DATA
else:
    LINCS_BC_DATA = BC_DATA / Path('data/data/raw/lincs')

lincs_metadata_path = LINCS_BC_DATA / "lincs_sig_info_new.csv"
LINCS_METADATA_PATH = LINCS_BC_DATA / "lincs_sig_info_new.csv"

BC_DISEASE_DATA = BC_DATA / "data" / "data"
BC_DISEASE_DIR = BC_DISEASE_DATA / DISEASE

############################################
# TSR directories
####################################
TSR_DIR=BASE_DIR / 'tsr'
TSR_OUT=TSR_DIR / 'output'
TSR_OUT_DISEASE= TSR_OUT / 'disease_signature' / DISEASE
TSR_OUT_DRUG=TSR_OUT / 'drug_signature'
TSR_OUT_DRUG_META=TSR_OUT_DRUG / 'LINCS_lorenzo' / 'metanalysis_mith3_drug_wise'
TSR_OUT_CSCORE=TSR_OUT / 'connectivity_score'

# Change here the data souce
DISEASE_SIGNATURE_SOURCE = "tsr"  # options: "tsr", "other"

############################################
# Connectivity score directories
############################################
CS_DIR=BASE_DIR / 'connectivity_score'
CS_IN_DRUG=CS_DIR / 'input' / 'drug_signature_2025' / cell_line_run_name #Path(cell_line+'_'+pert_time)
CS_IN_DISEASE=CS_DIR / 'input' / 'disease_signature_2025' / disease_run_name
CS_OUT=CS_DIR / 'output' / disease_run_name #Path(DISEASE+'_2025_'+pert_time)



# other outputs
IMG_DIR=BASE_DIR / 'imgs'
OOUT_DIR=BASE_DIR / 'other_outputs'


##########################################################################
# Drug signature calculations using LINCS1000 -
# Catalano dataset
##########################################################################

# LINCS_DIR=BASE_DIR / 'data' / 'LINCS-GSE92742'
# # LINCS1000 data filename (warning: big file):
# LINCS_file = LINCS_DIR / 'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
# #LINCS1000 metadata filename:
# inst_info_file = LINCS_DIR /  'GSE92742_Broad_LINCS_inst_info.txt'
# #LINCS1000 gene info filename:
# GENE_INFO_FILE = LINCS_DIR /  'GSE92742_Broad_LINCS_gene_info.txt'
# # landmark genes filename:
# BING_GENES = LINCS_DIR / "bing_gene_symbols.csv"  # optional



###############################################################################
# MITHrIL hyperparameters
###############################################################################
mith_batch_threads = 10
mith_organism = 'hsa'
mith_threads = 10

###############################################################################
#  Connectivity Score calculation hyperparameters and parameters
###############################################################################

# hyperparameters
cs_batch_threads = 1
cs_mith = 1 # [0,1] 1: calculate on MITHrIL data, 0: calculate on DEG data
cs_on_LM = 0 # possible values: [0,1] 0: calculate cs on all genes list 1: calculate on only landmark genes list
CS_METHOD ="bin_chen" #bin_chen_disease_sorted
CS_ON_PATHWAYS = True # Bool Default: False, calculate on signatures or pathways. Only works on mithril data

# Drug-collapse stage flag.
# False (BRCA/LINCS BinChen pipeline): the CS drug-input signatures are
#   individual LINCS1000 perturbations (one id per drug x time x dose, the 'ID'
#   column of the LINCS metadata, written as 'LINCS_id' in the cs_batch output).
#   pert_id (the drug name, with possible duplicates across timepoints/doses) is
#   mapped in AFTER the connectivity score via add_pert_id_to_cs, and drugs are
#   collapsed to one value per pert_id downstream.
# True (IPF/ALS_NYGC pipeline): drugs were already collapsed to one signature per drug via
#   a linear mixed model BEFORE MITHrIL and CS (see Catalano's thesis, koudjis2022).
#   The CS drug-input signature ids ARE pert_ids already, so the output id column
#   is written directly as 'pert_id' and add_pert_id_to_cs is skipped.
drug_collapsed_before_cs = True # Bool [True, False]

# CS functions


def validate_hyperparameters(cs_mith, cs_on_LM, CS_ON_PATHWAYS,
                             CS_METHOD=None, cs_threshold=None):
    """
    Raise ValueError if a hyperparameter combination is not supported.

    Used both at conf import time (against this file's defaults) and per-combo
    by the multi-conf wrappers (call_cs_batch_for_different_confs.py and
    hyperparameter_grid_search_bootstrap.py).

    The first three args are the pre-CS knobs and have always been required.
    CS_METHOD and cs_threshold are optional (default None = skip the check)
    so that legacy callers still work and the conf-import-time check below
    keeps using the original 3-arg form.
    """
    if cs_mith == 0 and CS_ON_PATHWAYS == 1:
        raise ValueError("Cannot load pathway signatures for DEG signatures."
                         " Set cs_mith to 1 and run MITHrIL propagations")
    if cs_on_LM == 1 and CS_ON_PATHWAYS == 1:
        raise ValueError("Cannot filter for landmark genes on pathways")

    # cs_threshold = -1.5 is calibrated for bin_chen_disease_sorted, whose
    # outputs are sums of normalized KS statistics. bin_chen returns RGES-style
    # values on a different scale, so a separate threshold must be calibrated
    # via RGES_study.py before bin_chen combos can be validated.
    # if (CS_METHOD == 'bin_chen' and cs_threshold is not None
    #         and float(cs_threshold) == -1.5):
    #     raise ValueError("cs_threshold=-1.5 is only meaningful for "
    #                      "CS_METHOD='bin_chen_disease_sorted'. For "
    #                      "CS_METHOD='bin_chen', calibrate a threshold first "
    #                      "(see RGES_study.py).")

# fail fast if the defaults in this conf are inconsistent
validate_hyperparameters(cs_mith, cs_on_LM, CS_ON_PATHWAYS)


def make_cs_filename(cs_mith, cs_on_LM, CS_ON_PATHWAYS, CS_METHOD, cs_out=None):
    """
    Build the timestamped output filename for a CS run.

    cs_out: optional Path. If None, falls back to module-level CS_OUT
        (i.e. the disease_run_name dictated by this conf.py). The wrapper
        passes its own cs_out so non-default landmark_disease lands in the
        right directory.
    """
    if cs_out is None:
        cs_out = CS_OUT
    now = datetime.now()
    datetime_string = now.strftime("%d_%m_%Y_%H_%M")
    mith = "_mith" if cs_mith == 1 else "_DEG"
    LM = "_LM" if cs_on_LM == 1 else ""
    pathways = "_pw" if CS_ON_PATHWAYS == 1 else ""
    score = "_"+CS_METHOD
    cs_filename = datetime_string + mith + LM + pathways + score + "_connectivity_score.tsv"

    connectivity_dataset_filename = cs_out / cs_filename
    cs_id = Path(cs_filename).stem
    return connectivity_dataset_filename, cs_id

# log filename
cs_log_filename = Path(LOGS_DIR) / 'cs_runs.tsv'


###############################################################################
# BinChen 2017 sRGES replication parameters
# (validations/binchen2017/test_sRGES_replication.py)
###############################################################################

# cell line per disease
DEFAULT_CELL_LINE = {
    "LIHC": "HEPG2",
    "BRCA": "MCF7",
    "COAD": "HT29",
}

# sRGES algorithmic mode: 'binchen_r' faithfully reproduces the R source;
# 'corrected' fixes the integer-index quirk (see src/sRGES.py docstring).
srges_diff_mode = "binchen_r"

# Apply per-cell-line cor*RGES weighting (Bin Chen sRGES_all_cmpds.R).
srges_use_cor = False

# Default output directory for binchen2017 validation reports.
BINCHEN_VALIDATION_OUT_DIR = OOUT_DIR / "binchen2017_validation"


###############################################################################
# Chembl IC50 validation parameters
###############################################################################

# FIles and directories
VAL_DIR = BASE_DIR / 'validations' 

# Chembl directory
CHEMBL_BASE_DIR = VAL_DIR / 'chembl'

CHEMBL_INPUT_DATA_DIR = CHEMBL_BASE_DIR / 'chembl_input'
# cell_lines_chembl = ["MCF7", "HepG2", "HT29"]
ic50_file = DATA_DIR/'BinChen2017'/'SD8.xlsx'

chembl_val_log_filename = LOGS_DIR / "BinChen2017_chembl_validation_IC50_correlation_runs.tsv"

# Hyperparameters:

IC_50_binchen_SD5= True # to implement the False: load freshly downloaded chembl  

CS_DRUG_COLLAPSE_METHOD= 'srges'#'best' # None #
IC50_DRUG_COLLAPSE_METHOD='median' # 

# IC50 data filtering
IC50_ONLY=True

# Thresholds
CS_TH = -1.5  # Connectivity score threshold
IC50_EFF_TH = 10.0 # IC50 effective threshold


# ------------------------------------------------------------------
# Optional: manual selection of an existing CS run for downstream steps
# Leave both as None for automatic resolution from cs_runs.tsv
# ------------------------------------------------------------------
selected_cs_run_id = None


####################
# Utility functions
####################

# load map of gene symols : gene ids
# more symbols may map to the same gene id

# using signorini's mapping (see network_signing algorithm)
alias_2geneid_filename=DICT_DIR / 'alias_2geneid.pkl'
with open(alias_2geneid_filename, 'rb') as f:
    alias_2geneid = pickle.load(f)

def map_name_to_id(genename):
    if genename in alias_2geneid.keys():
        return alias_2geneid[genename]
    return genename

# # using alaimo's mapping
# symbol_to_id_filename=BASE_DIR / 'other_data' / 'symbol_2geneid.pkl'
# with open(symbol_to_id_filename, 'rb') as f:
#     symbol_2geneid = pickle.load(f)

# def map_name_to_id(genename):
#     if genename in symbol_2geneid.keys():
#         return symbol_2geneid[genename]
#     return genename
