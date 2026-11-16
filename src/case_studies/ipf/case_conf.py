# -*- coding: utf-8 -*-
"""
Specification of the ipf case study: the standard paths from
conf.case_study_paths (slug = name of this folder) and the file names of the
GSE92592 data (IPF and control lung tissue RNA-seq).

The raw counts behind the ipf disease signature are in conf.TSR_DIR / 'data'
/ 'IPF-GSE92592'; the NCBI TPM table and the series matrix are in the
case-study data folder (DATA_DIR / slug).

The analysis fields at the bottom (SAMPLES, FIGURES and the validation
declarations) are read by case_study_pipeline.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from conf import case_study_paths, TSR_DIR

P = case_study_paths(Path(__file__).resolve().parent.name)

NCBI_COUNTS = TSR_DIR / 'data' / 'IPF-GSE92592' / 'GSE92592_raw_counts_GRCh38.p13_NCBI.tsv'
SERIES_MATRIX = P.data_dir / 'GSE92592_series_matrix.txt.gz'
NCBI_TPM = P.data_dir / 'GSE92592_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz'

# sample order and group coding of the disease signature (0 = IPF, 1 = control)
SAMPLES_GROUPS_MAP = '000000000000000000001111111111111111111'

# ---------------------------------------------------------------------------
# Analysis specification read by case_study_pipeline.py (fields not declared
# here = analysis not run for this study).

# step 'weights': disease samples; the first group is the case group
SAMPLES = {
    'format': 'series_matrix',
    'file': SERIES_MATRIX,
    'characteristic': 'disease state',   # !Sample_characteristics_ch1 field of the group
    'control_value': 'Control',          # its value for controls; any other value = case
    'groups': {'ipf': '0', 'control': '1'},   # group: group_code
    'check_group_coding': SAMPLES_GROUPS_MAP,
    'check_sample_order': NCBI_COUNTS,
}

# step 'validation': known drugs (name: class)
KNOWN_DRUGS = {
    'pirfenidone': 'approved, disease-modifying',
    'nintedanib': 'approved, disease-modifying',
}

# step 'validation': known-answer checks of the per-instance CS of every configuration
CS_CHECKS = {
    'cell_line': 'NEU',             # cell line of the recomputation check
}

# figure text
FIGURES = {
    'dataset': 'GSE92592',
    'disease': 'IPF',
    'group_labels': {'ipf': 'IPF'},
    'case_axis': 'IPF',
    'control_axis': 'controls',
    'panel_c_title': 'IPF vs control samples',
    'known_drugs_label': 'IPF',
}
