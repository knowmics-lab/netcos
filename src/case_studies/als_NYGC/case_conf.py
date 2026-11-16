# -*- coding: utf-8 -*-
"""
Specification of the als_NYGC case study: the standard paths from
conf.case_study_paths (slug = name of this folder, = conf.DISEASE id of the
als_NYGC connectivity-score runs) and the file names of GSE153960 (NYGC ALS
Consortium / Target ALS postmortem CNS RNA-seq).

The NCBI TPM matrix is in the slug data folder; the raw counts and the sample
metadata behind the als_NYGC disease signature are in the TSR data folder.

The analysis fields at the bottom (SAMPLES, FIGURES and the validation
declarations) are read by case_study_pipeline.py.
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

# ---------------------------------------------------------------------------
# Analysis specification read by case_study_pipeline.py (fields not declared
# here = analysis not run for this study).

# step 'weights': disease samples; the first group is the case group
SAMPLES = {
    'format': 'csv',
    'file': METADATA,
    'read_csv': {'sep': ';', 'encoding': 'utf-8-sig'},
    'sample_column': 'accession',
    'group_column': 'tissue_status',
    'groups': {'als': CASE_STATUS, 'control': CONTROL_STATUS},   # group: value of group_column
    'columns': ['GSM', 'tissue_status', 'group', 'tissue', 'region_class', 'subject_id'],
    'subset_of_tpm': True,          # the TPM matrix also holds samples of other groups
    'strata': {                     # CNS region class of a sample, from its tissue
        'name': 'region_class',
        'from_column': 'tissue',
        'prefixes': {'cortex': 'cortex', 'spinal cord': 'spinal_cord'},
        'classes': ['cortex', 'spinal_cord', 'cerebellum', 'hippocampus'],
        'label': 'CNS region class',
        'short_label': 'region class',
    },
}

# step 'validation': known drugs (name: class) and the known-drug report
KNOWN_DRUGS = {
    'riluzole': 'approved, disease-modifying',
    'edaravone': 'approved, disease-modifying',
    'dextromethorphan': 'approved, symptomatic (Nuedexta component)',
    'quinidine': 'approved, symptomatic (Nuedexta component)',
    'phenylbutyrate': 'failed phase 3 (AMX0035 component)',
}
KNOWN_DRUGS_REPORT = 'percentiles'

# step 'validation': known-answer checks of the per-instance DEG CS
CS_CHECKS = {
    'disease_signature': DISEASE_SIGNATURE,
    'cell_line': 'NEU',             # cell line of the recomputation check
    'n_check': 20,                  # instances recomputed
    'n_jobs': 12,                   # parallel jobs of the position-convention check
}

# figure text and style
FIGURES = {
    'style': 'lfs',
    'rc': {'savefig.dpi': 130, 'font.size': 10, 'legend.shadow': False},
    'dataset': 'GSE153960',
    'disease': 'ALS CNS',
    'case_axis': 'als samples',
    'control_axis': 'control samples',
    'panel_c_title': 'als vs control samples',
    'panel_a_title': 'LINCS lines with CCLE expression',
    'panel_a_legend': 'upper right',
    'known_drugs_label': 'ALS',
}
