#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-instance drug inputs for every LINCS cell line, from the Chen 2017 L1000
level-5 landmark signatures: direct connectivity-score inputs and MITHrIL
input matrices.

Design
------
- Signatures: lincs_signatures_cmpd_landmark.RData (rows = 978 landmark NCBI
  GeneIDs, columns = LINCS instance ids); instances and their cell line from
  lincs_sig_info_new.csv (all rows of a cell line, no further filter).
- One worker per cell line (multiprocessing; the signature matrix is shared
  by fork). Every cell line gets the run name '<cell_id>_LM' (landmark genes
  only).
- Existing files are never overwritten: a CS pickle already present is
  skipped, and the MITHrIL matrix of a line is written only if absent.
- The run registry row of each line is appended by the parent process after
  all workers finish.

Inputs
------
- conf.LINCS_BC_DATA / 'lincs_signatures_cmpd_landmark.RData'
- conf.LINCS_METADATA_PATH (lincs_sig_info_new.csv: id, cell_id, pert_iname, ...)

Outputs
-------
- conf.CS_DIR/input/drug_signature_2025/<cell_id>_LM/<id>_signature_gene_id.pkl
    one DataFrame per instance: signature_id (str), gene_id (str),
    DE_log2_FC (level-5 value), adj.p.value (= 1.0), 978 rows.
- conf.MITH_IN_DRUG/LINCS_<cell_id>_LM.mi
    genes x instances, tab-separated; header = instance ids (no leading
    cell), then one row per gene: gene_id followed by the values.
- conf.LOGS_DIR/lincs_drug_preprocessing_runs.tsv
    one row per cell line processed.

Usage
-----
    cd src && PYTHONPATH=. python preprocess_binchen2017_lincs_drug_inputs.py
Set ONLY_CELL_LINES to a list of cell_id to process a subset.
"""

import os
import socket
from datetime import datetime
from multiprocessing import Pool

import pandas as pd
import pyreadr

from conf import LINCS_BC_DATA, LINCS_METADATA_PATH, CS_DIR, MITH_IN_DRUG, LOGS_DIR
from logger import append_run_metadata

RUN_SUFFIX = '_LM'
N_WORKERS = 24
ONLY_CELL_LINES = None  # e.g. ['HS27A']

LINCS_RDATA = LINCS_BC_DATA / 'lincs_signatures_cmpd_landmark.RData'
CS_DRUG_ROOT = CS_DIR / 'input' / 'drug_signature_2025'
REGISTRY = LOGS_DIR / 'lincs_drug_preprocessing_runs.tsv'

FC = None  # genes x instances; set by load_inputs() and inherited by the workers
MD = None


def load_inputs():
    global FC, MD
    FC = pyreadr.read_r(LINCS_RDATA)['lincs_signatures']
    FC.index = FC.index.astype(str)
    FC.columns = FC.columns.astype(str)
    MD = pd.read_csv(LINCS_METADATA_PATH)
    MD['id'] = MD['id'].astype(str)


def process_cell_line(cell_id, cs_root=CS_DRUG_ROOT, mith_dir=MITH_IN_DRUG):
    run_name = cell_id + RUN_SUFFIX
    ids = [i for i in MD.loc[MD['cell_id'] == cell_id, 'id'] if i in FC.columns]
    cs_dir = cs_root / run_name
    cs_dir.mkdir(parents=True, exist_ok=True)
    n_written = 0
    for i in ids:
        out = cs_dir / f'{i}_signature_gene_id.pkl'
        if out.exists():
            continue
        df = pd.DataFrame({'signature_id': i, 'gene_id': FC.index,
                           'DE_log2_FC': pd.to_numeric(FC[i], errors='coerce'),
                           'adj.p.value': 1.0}).dropna(subset=['DE_log2_FC'])
        df.reset_index(drop=True).to_pickle(out)
        n_written += 1
    mi = mith_dir / f'LINCS_{run_name}.mi'
    mi_written = False
    if not mi.exists():
        tmp = mi.with_suffix('.mi.tmp')
        with open(tmp, 'w') as f:
            f.write('\t'.join(ids) + '\n')
            FC[ids].to_csv(f, sep='\t', header=False)
        os.replace(tmp, mi)
        mi_written = True
    return {'cell_line': cell_id, 'cell_line_run_name': run_name,
            'n_instances': len(ids),
            'n_unique_compounds': MD.loc[MD['id'].isin(ids), 'pert_iname'].nunique(),
            'n_cs_files_written': n_written, 'mith_input_written': int(mi_written),
            'cs_input_dir': str(cs_dir), 'mith_input_file': str(mi)}


if __name__ == '__main__':
    load_inputs()
    print(f'signatures: {FC.shape[0]} genes x {FC.shape[1]} instances')
    print(f'metadata: {len(MD)} rows, {MD["cell_id"].nunique()} cell lines; '
          f'ids missing from the signatures: {(~MD["id"].isin(FC.columns)).sum()}')
    cells = ONLY_CELL_LINES or sorted(MD['cell_id'].unique())
    MITH_IN_DRUG.mkdir(parents=True, exist_ok=True)
    with Pool(min(N_WORKERS, len(cells))) as pool:
        results = pool.map(process_cell_line, cells)
    timestamp = datetime.now().strftime('%d_%m_%Y_%H_%M')
    for r in results:
        print(r['cell_line'], r['n_instances'], 'instances,', r['n_cs_files_written'], 'CS files written,',
              'mi written' if r['mith_input_written'] else 'mi present')
        row = {'preprocessing_run_id': f'{timestamp}__{r["cell_line_run_name"]}', 'timestamp': timestamp,
               'hostname': socket.gethostname(), 'lincs_signatures_source_file': str(LINCS_RDATA),
               'lincs_metadata_source_file': str(LINCS_METADATA_PATH), **r}
        append_run_metadata(REGISTRY, row)
