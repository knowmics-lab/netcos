# -*- coding: utf-8 -*-
"""
Locations of the asthma_GSE152004 case study: the standard paths from
conf.case_study_paths (slug = name of this folder) and the file names of the
GSE152004 download.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from conf import case_study_paths

P = case_study_paths(Path(__file__).resolve().parent.name)

NCBI_COUNTS = P.data_dir / 'GSE152004_raw_counts_GRCh38.p13_NCBI.tsv.gz'
SERIES_MATRIX = P.data_dir / 'GSE152004_series_matrix.txt.gz'
