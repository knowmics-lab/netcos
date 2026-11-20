# -*- coding: utf-8 -*-
"""
Specification of the asthma_GSE152004 case study: the standard paths from
conf.case_study_paths (slug = name of this folder) and the file names of the
GSE152004 download.

The analysis fields at the bottom (SAMPLES, FIGURES and the validation
declarations) are read by case_study_pipeline.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from conf import case_study_paths

P = case_study_paths(Path(__file__).resolve().parent.name)

NCBI_COUNTS = P.data_dir / 'GSE152004_raw_counts_GRCh38.p13_NCBI.tsv.gz'
SERIES_MATRIX = P.data_dir / 'GSE152004_series_matrix.txt.gz'
NCBI_TPM = P.data_dir / 'GSE152004_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz'

# ---------------------------------------------------------------------------
# Analysis specification read by case_study_pipeline.py (fields not declared
# here = analysis not run for this study).

# step 'weights': disease samples; the first group is the case group
SAMPLES = {
    'format': 'samples_table',      # P.oout_dir/<slug>_samples.tsv of 00_inspect_geo_download.py
    'group_column': 'group_code',
    'groups': {'asthmatic': '0', 'control': '1'},   # group: value of group_column
}

# step 'validation': AUROC of approved vs negative repoDB drugs of the indication
REPODB = {
    'cui': 'C0004096',              # repoDB ind_id of asthma
    'negative_classes': ['efficacy', 'safety', 'efficacy_and_safety'],
    'negative_sets': {              # name: column whose value must be in negative_classes
        'stop_efficacy_safety': 'pair_stop_class',
        'stop_efficacy_safety_enrolled': 'clinical_stop_enrolled',
    },
    'name_alias': {'beclometasone': 'beclomethasonedipropionate'},   # repoDB name key -> L1000 name key
}

# step 'validation': recomputation checks of the per-instance CS of every configuration
CS_CHECKS = {
    'cell_line': 'NEU',             # cell line of the recomputation check
}

# figure text
FIGURES = {
    'dataset': 'GSE152004',
    'disease': 'asthma',
    'case_axis': 'asthmatics',
    'control_axis': 'healthy controls',
    'panel_c_title': 'asthmatic vs healthy-control samples',
}
