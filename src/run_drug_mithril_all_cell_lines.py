#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MITHrIL on the per-instance LINCS drug signatures of every cell line, followed
by the conversion of its outputs to connectivity-score inputs (gene level and
pathway level).

Design
------
- Cell lines are processed one at a time, from the smallest to the largest
  number of instances; each line is one mithril-batch run with N_THREADS
  threads and the parameters of step6_drug_run_MITHrIL.run_mithril_batch
  (-reactome, -customize-pathway-matrix, -seed 1234, fast-cpu factories,
  organism conf.mith_organism).
- Input: conf.MITH_IN_DRUG/LINCS_<cell_id>_LM.mi. If part of a line's
  instances already have an output, a matrix with the remaining instances only
  (LINCS_<cell_id>_LM.remaining.mi) is run instead, so an interrupted run
  resumes where it stopped.
- After MITHrIL, every <id>.perturbations.txt is converted with
  step7_drug_mith_out_to_cs_in.convert_directory (-> <id>.pkl) and every
  <id>.output.txt with step2b_7b_drug_and_disease_mith_out_to_cs_in.
  convert_one_file (-> <id>_pathways.tsv). Existing outputs are not rewritten.
- A line is skipped when all its instances have both CS inputs. MITHrIL
  outputs count as present also when zstd-compressed (<file>.zst, see
  compress_mithril_drug_outputs.py); step7 reads them directly.

Inputs
------
- conf.MITH_IN_DRUG/LINCS_<cell_id>_LM.mi (from preprocess_binchen2017_lincs_drug_inputs.py)
- conf.LINCS_METADATA_PATH (instances per cell line)
- conf.MITH_APP, conf.mith_organism

Outputs
-------
- conf.MITH_OUT/drug_signature_2025/<cell_id>_LM/<id>.perturbations.txt, <id>.output.txt
- conf.CS_DIR/input/drug_signature_2025/<cell_id>_LM/<id>.pkl, <id>_pathways.tsv
- conf.LOGS_DIR/mithril_drug_runs.tsv: one row per mithril-batch run
  (wall-clock elapsed_sec, return code, instances run, outputs found).

Usage
-----
    cd src && PYTHONPATH=. python run_drug_mithril_all_cell_lines.py
Set ONLY_CELL_LINES to a list of cell_id to process a subset.
"""

import socket
import subprocess
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

from loader import read_lincs_metadata

from conf import (CS_DIR, LINCS_METADATA_PATH, LOGS_DIR, MITH_APP, MITH_IN_DRUG, MITH_OUT,
                  mith_organism)
from logger import append_run_metadata
from step7_drug_mith_out_to_cs_in import convert_directory
from step2b_7b_drug_and_disease_mith_out_to_cs_in import convert_one_file

RUN_SUFFIX = '_LM'
N_THREADS = 28
ONLY_CELL_LINES = None

MITH_OUT_ROOT = MITH_OUT / 'drug_signature_2025'
CS_DRUG_ROOT = CS_DIR / 'input' / 'drug_signature_2025'
REGISTRY = LOGS_DIR / 'mithril_drug_runs.tsv'


def mithril_command(in_file, out_dir, n_threads):
    return ['java', '-jar', str(MITH_APP), 'mithril-batch', '-reactome', '-p', '-customize-pathway-matrix',
            '-seed', '1234', '-inversion-factory', 'fast-cpu', '-multiplication-factory', 'fast-cpu',
            '-t', str(n_threads), '-organism', mith_organism, '-i', str(in_file), '-o', str(out_dir),
            '-verbose']


def done_ids(out_dir, suffix):
    """Ids with a file '<id><suffix>' (or its zstd-compressed '<id><suffix>.zst') in out_dir."""
    if not out_dir.exists():
        return set()
    ids = set()
    for sfx in (suffix, suffix + '.zst'):
        ids |= {p.name[:-len(sfx)] for p in out_dir.glob('*' + sfx)
                if not p.name.endswith('_signature_gene_id.pkl')}
    return ids


def run_cell_line(cell_id, ids):
    run_name = cell_id + RUN_SUFFIX
    mi = MITH_IN_DRUG / f'LINCS_{run_name}.mi'
    out_dir = MITH_OUT_ROOT / run_name
    cs_dir = CS_DRUG_ROOT / run_name
    out_dir.mkdir(parents=True, exist_ok=True)

    have = done_ids(out_dir, '.perturbations.txt') & done_ids(out_dir, '.output.txt')
    todo = [i for i in ids if i not in have]
    if todo:
        in_file = mi
        if len(todo) < len(ids):
            fc = pd.read_csv(mi, sep='\t')  # header has one field fewer: first column becomes the index
            in_file = MITH_IN_DRUG / f'LINCS_{run_name}.remaining.mi'
            with open(in_file, 'w') as f:
                f.write('\t'.join(todo) + '\n')
                fc[todo].to_csv(f, sep='\t', header=False)
        cmd = mithril_command(in_file, out_dir, N_THREADS)
        print(f'[{datetime.now():%d/%m %H:%M}] {run_name}: MITHrIL on {len(todo)} of {len(ids)} instances', flush=True)
        start = time.time()
        rc = subprocess.run(cmd).returncode
        elapsed = time.time() - start
        n_out = len(done_ids(out_dir, '.perturbations.txt') & set(ids))
        ts = datetime.now().strftime('%d_%m_%Y_%H_%M_%S')
        append_run_metadata(REGISTRY, {
            'mithril_run_id': f'{ts}__{run_name}', 'timestamp': ts, 'hostname': socket.gethostname(),
            'step_name': 'run_drug_mithril_all_cell_lines', 'run_scope': 'drug_batch',
            'drug_run_id': run_name, 'cell_line': cell_id,
            'filter_drug_signature_by_landmark_genes_pre_mith': 1,
            'mith_input_file': str(in_file), 'mith_output_dir': str(out_dir),
            'mith_app': str(MITH_APP), 'mith_organism': mith_organism, 'mith_threads': N_THREADS,
            'return_code': rc, 'elapsed_sec': round(elapsed, 1), 'n_instances_run': len(todo),
            'n_drugs_output': n_out, 'n_instances_cell_line': len(ids)})
        print(f'[{datetime.now():%d/%m %H:%M}] {run_name}: return code {rc}, {elapsed / 3600:.2f} h, '
              f'{n_out}/{len(ids)} outputs', flush=True)
        if rc != 0:
            return

    cs_dir.mkdir(parents=True, exist_ok=True)
    if done_ids(out_dir, '.perturbations.txt') - done_ids(cs_dir, '.pkl'):
        convert_directory(out_dir, cs_dir, drop_special_gene_names=False)
    for i in done_ids(out_dir, '.output.txt') - done_ids(cs_dir, '_pathways.tsv'):
        convert_one_file(out_dir / f'{i}.output.txt', cs_dir / f'{i}_pathways.tsv')
    print(f'[{datetime.now():%d/%m %H:%M}] {run_name}: CS inputs '
          f'{len(done_ids(cs_dir, ".pkl"))} gene-level, {len(done_ids(cs_dir, "_pathways.tsv"))} pathway', flush=True)


if __name__ == '__main__':
    md = read_lincs_metadata(LINCS_METADATA_PATH)
    md['id'] = md['id'].astype(str)
    sizes = md.groupby('cell_id').size().sort_values()
    cells = ONLY_CELL_LINES or sizes.index.tolist()
    for c in cells:
        ids = md.loc[md['cell_id'] == c, 'id'].tolist()
        if not (MITH_IN_DRUG / f'LINCS_{c}{RUN_SUFFIX}.mi').exists():
            print(f'{c}: no MITHrIL input, skipped', flush=True)
            continue
        hdr = open(MITH_IN_DRUG / f'LINCS_{c}{RUN_SUFFIX}.mi').readline().rstrip('\n').split('\t')
        run_cell_line(c, hdr)
