#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-drug scores of GSE152004 from the per-instance connectivity scores of all
LINCS cell lines (Chen 2017 / OCTAD route: drugs not collapsed before CS).

Design
------
- Per-instance CS: for every configuration in CONFIGS and every LINCS cell
  line, the most recent cs_batch run of conf.DISEASE on drug input
  '<cell_id>_LM', resolved with loader.load_cs_run; the cell lines with no run
  for a configuration are listed and left out.
- Instance metadata (pert_id, pert_iname, cell_id, pert_dose, pert_time) from
  conf.LINCS_METADATA_PATH, joined on the LINCS instance id.
- Three per-drug collapses (drug = pert_id, the L1000 compound id):
  * 'sRGES_weighted': sRGES.compute_sRGES with the cell-line weights of
    02_cell_line_weights.py (cor = median correlation with the asthmatic
    samples); instances of lines without a weight are left out.
  * 'sRGES_unweighted': sRGES.compute_sRGES with every weight = 1, all lines.
  * 'best': the most negative CS among the drug's instances with
    pert_dose > 0 and pert_time in {6, 24} h (the instances sRGES keeps).
  sRGES uses the Chen 2017 dose/time correction (srges_diff_mode 'binchen_r').
- Sign convention: more negative = stronger reversal of the disease signature.

Inputs
------
- conf.cs_log_filename and the CS files it points to
- conf.LINCS_METADATA_PATH
- P.oout_dir/<slug>_cell_line_weights.tsv

Outputs
-------
- P.oout_dir/<slug>_srges_per_drug.tsv
    one row per (config, collapse, pert_id): pert_iname, score, n_instances,
    n_cell_lines.
- P.oout_dir/<slug>_srges_cs_runs_used.tsv
    the cs_run_id used for every (config, cell line).

Usage
-----
    cd src && PYTHONPATH=. python case_studies/asthma_GSE152004/03_srges_per_drug.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from case_conf import P
import conf
from loader import load_cs_run
from sRGES import compute_sRGES

CONFIGS = {  # label: (mith, cs_on_LM, CS_ON_PATHWAYS)
    'DEG': (0, 0, 0),
    'mith': (1, 0, 0),
    'mith_LM': (1, 1, 0),
    'mith_pw': (1, 0, 1),
}
RUN_SUFFIX = '_LM'


def load_instances(label, cells):
    m, lm, pw = CONFIGS[label]
    parts, used, missing = [], [], []
    for c in cells:
        try:
            df, path, run_id = load_cs_run(disease=conf.DISEASE, landmark_disease=False,
                                           drug_run_id=c + RUN_SUFFIX, mith=m, cs_on_LM=lm,
                                           CS_ON_PATHWAYS=pw, CS_METHOD=conf.CS_METHOD, verbose=False)
        except FileNotFoundError:
            missing.append(c)
            continue
        parts.append(df[['LINCS_id', 'connectivity_score']])
        used.append({'config': label, 'cell_line': c, 'cs_run_id': run_id, 'n_rows': len(df)})
    inst = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=['LINCS_id', 'connectivity_score'])
    return inst, used, missing


def collapse(inst, md, weights):
    x = inst.assign(LINCS_id=inst['LINCS_id'].astype(str)).merge(md, left_on='LINCS_id', right_on='id')
    x['pert_dose'] = pd.to_numeric(x['pert_dose'], errors='coerce')
    x['pert_time'] = pd.to_numeric(x['pert_time'], errors='coerce')
    kept = x[(x['pert_dose'] > 0) & x['pert_time'].isin([6, 24])]
    names = md.drop_duplicates('pert_id').set_index('pert_id')['pert_iname']
    out = []
    for name, w in [('sRGES_weighted', weights), ('sRGES_unweighted', None)]:
        s = compute_sRGES(x, score_col='connectivity_score', drug_col='pert_id', cell_col='cell_id',
                          cell_line_weights=w, srges_diff_mode='binchen_r')
        lines = (kept[kept['cell_id'].isin(weights)] if w is not None else kept).groupby('pert_id')['cell_id'].nunique()
        out.append(pd.DataFrame({'collapse': name, 'pert_id': s['pert_id'], 'score': s['sRGES'],
                                 'n_instances': s['n'], 'n_cell_lines': s['pert_id'].map(lines)}))
    g = kept.groupby('pert_id')
    out.append(pd.DataFrame({'collapse': 'best', 'pert_id': g['connectivity_score'].min().index,
                             'score': g['connectivity_score'].min().values,
                             'n_instances': g.size().values, 'n_cell_lines': g['cell_id'].nunique().values}))
    res = pd.concat(out, ignore_index=True)
    res.insert(2, 'pert_iname', res['pert_id'].map(names))
    return res


if __name__ == '__main__':
    md = pd.read_csv(conf.LINCS_METADATA_PATH, dtype={'id': str})
    md = md[['id', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose', 'pert_time']]
    cells = sorted(md['cell_id'].unique())
    w = pd.read_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t')
    weights = dict(zip(w['lincs_cell_id'], w['cor']))
    print(f'{len(cells)} LINCS cell lines; {len(weights)} with a weight')

    tabs, used_all = [], []
    for label in CONFIGS:
        inst, used, missing = load_instances(label, cells)
        used_all += used
        print(f'{label}: {len(used)} cell lines, {len(inst)} instances; '
              f'no run for {len(missing)} lines' + (f': {", ".join(missing)}' if missing else ''))
        if inst.empty:
            continue
        r = collapse(inst, md, weights)
        r.insert(0, 'config', label)
        r['n_cell_lines_with_run'] = len(used)
        tabs.append(r)
        print(r.groupby('collapse').agg(n_drugs=('pert_id', 'size'), median_score=('score', 'median')).round(4)
              .to_string())

    P.oout_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(tabs, ignore_index=True).to_csv(P.oout_dir / f'{P.slug}_srges_per_drug.tsv', sep='\t', index=False)
    pd.DataFrame(used_all).to_csv(P.oout_dir / f'{P.slug}_srges_cs_runs_used.tsv', sep='\t', index=False)
