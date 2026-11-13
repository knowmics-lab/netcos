#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-instance connectivity scores of one disease against the LINCS drug
instances of every cell line, one cs_batch run per (cell line, configuration).

Design
------
- Drug inputs: conf.CS_DIR/input/drug_signature_2025/<cell_id>_LM/ (one file
  per LINCS instance, drugs not collapsed before CS).
- Disease: conf.DISEASE, disease input directory as in cs_batch.
- Configurations: CONFIGS, each a (cs_mith, cs_on_LM, CS_ON_PATHWAYS) triple,
  scored with conf.CS_METHOD. A MITHrIL configuration is run on a cell line
  only when every instance of the line has its MITHrIL-derived input.
- A (cell line, configuration) already present in conf.cs_log_filename for
  this disease is skipped. Runs start in different minutes, so the
  minute-stamped CS file names of cs_batch never collide.

Inputs
------
- conf.DISEASE, conf.CS_METHOD, conf.LINCS_METADATA_PATH, conf.cs_log_filename
- drug and disease CS inputs (see above)

Outputs
-------
- conf.CS_DIR/output/<DISEASE>/<timestamp>_<config>_<CS_METHOD>_connectivity_score.tsv
  (written by cs_batch.run_cs_batch_for_conf, one per run) and one row per run
  in conf.cs_log_filename (drug_run_id = <cell_id>_LM).

Usage
-----
    cd src && PYTHONPATH=. python run_cs_all_cell_lines.py
"""

import time
from datetime import datetime

import pandas as pd

import conf
from cs_batch import run_cs_batch_for_conf
from loader import read_lincs_metadata

RUN_SUFFIX = '_LM'
N_JOBS = 6
ONLY_CELL_LINES = None
CONFIGS = {  # label: (cs_mith, cs_on_LM, CS_ON_PATHWAYS)
    'DEG': (0, 0, False),
    'mith': (1, 0, False),
    'mith_LM': (1, 1, False),
    'mith_pw': (1, 0, True),
}
RUN_CONFIGS = ['DEG', 'mith', 'mith_LM', 'mith_pw']


def mith_inputs_complete(cs_dir, ids, pathways):
    suffix = '_pathways.tsv' if pathways else '.pkl'
    present = {p.name[:-len(suffix)] for p in cs_dir.glob('*' + suffix)}
    return all(i in present for i in ids)


def already_run(log, run_name, cfg):
    if log is None:
        return False
    m, lm, pw = cfg
    sel = log[(log['disease_run_id'] == conf.DISEASE) & (log['drug_run_id'] == run_name)
              & (log['mith'] == m) & (log['cs_on_LM'] == lm)
              & (log['CS_ON_PATHWAYS'].astype(str).isin([str(pw), str(int(pw))]))
              & (log['CS_METHOD'] == conf.CS_METHOD) & (log['drug_collapsed_before_cs'] == 0)]
    return len(sel) > 0


if __name__ == '__main__':
    md = read_lincs_metadata()
    md['id'] = md['id'].astype(str)
    cells = ONLY_CELL_LINES or md.groupby('cell_id').size().sort_values().index.tolist()
    last_minute = None
    for label in RUN_CONFIGS:
        cfg = CONFIGS[label]
        for c in cells:
            run_name = c + RUN_SUFFIX
            cs_dir = conf.CS_DIR / 'input' / 'drug_signature_2025' / run_name
            ids = md.loc[md['cell_id'] == c, 'id'].tolist()
            log = pd.read_csv(conf.cs_log_filename, sep='\t') if conf.cs_log_filename.exists() else None
            if already_run(log, run_name, cfg):
                continue
            if cfg[0] == 1 and not mith_inputs_complete(cs_dir, ids, cfg[2]):
                print(f'{label} {run_name}: MITHrIL inputs incomplete, skipped', flush=True)
                continue
            while datetime.now().strftime('%d%m%Y%H%M') == last_minute:
                time.sleep(5)
            last_minute = datetime.now().strftime('%d%m%Y%H%M')
            print(f'[{datetime.now():%d/%m %H:%M}] {label} {run_name} ({len(ids)} instances)', flush=True)
            run_cs_batch_for_conf(cell_line=c, landmark_drug=True, landmark_disease=False,
                                  CS_METHOD=conf.CS_METHOD, cs_mith=cfg[0], cs_on_LM=cfg[1],
                                  CS_ON_PATHWAYS=cfg[2], cs_batch_threads=min(N_JOBS, len(ids)),
                                  DISEASE=conf.DISEASE, drug_collapsed=False)
