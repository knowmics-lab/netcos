# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:11:23 2026

@author: los4
"""

from cs_batch import run_cs_batch_for_conf
from conf import CS_IN_DISEASE, CS_DIR, DISEASE
CS_IN_DRUG = CS_DIR/'input'/'drug_signature_2025'/'all_LMM_metanalysis'

run_cs_batch_for_conf(CS_IN_DISEASE=CS_IN_DISEASE, CS_IN_DRUG=CS_IN_DRUG, DISEASE=DISEASE)
