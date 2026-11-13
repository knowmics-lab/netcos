# -*- coding: utf-8 -*-
"""
Locations of the ipf case study: the standard paths from
conf.case_study_paths (slug = name of this folder) and the file names of the
GSE92592 data (IPF and control lung tissue RNA-seq).

The raw counts behind the ipf disease signature are in conf.TSR_DIR / 'data'
/ 'IPF-GSE92592'; the NCBI TPM table and the series matrix are in the
case-study data folder (DATA_DIR / slug).
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
