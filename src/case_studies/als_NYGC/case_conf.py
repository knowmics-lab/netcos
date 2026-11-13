# -*- coding: utf-8 -*-
"""
Locations of the als_NYGC case study: the standard paths from
conf.case_study_paths (slug = name of this folder, = conf.DISEASE id of the
als_NYGC connectivity-score runs) and the file names of GSE153960 (NYGC ALS
Consortium / Target ALS postmortem CNS RNA-seq).

The NCBI TPM matrix is in the slug data folder; the raw counts and the sample
metadata behind the als_NYGC disease signature are in the TSR data folder.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from conf import case_study_paths, TSR_DIR, TSR_OUT

P = case_study_paths(Path(__file__).resolve().parent.name)

NCBI_TPM = P.data_dir / 'GSE153960_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz'
TSR_DATA_DIR = TSR_DIR / 'data' / 'ALS-NYGC-GSE153960'
METADATA = TSR_DATA_DIR / 'gse153960_metadata.csv'
NCBI_COUNTS = TSR_DATA_DIR / 'GSE153960_raw_counts_GRCh38.p13_NCBI.tsv.gz'

# DEG disease signature of the als_NYGC connectivity-score runs (TSR / LMM output)
DISEASE_SIGNATURE = TSR_OUT / 'disease_signature' / P.slug / f'{P.slug}_signature_gene_id.csv'

# sample groups of the als_NYGC disease signature (metadata column tissue_status)
CASE_STATUS = 'ALS Spectrum MND'
CONTROL_STATUS = 'Non-Neurological Control'
