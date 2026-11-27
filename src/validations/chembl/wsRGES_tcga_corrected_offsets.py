#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-drug scores of BRCA, LIHC and COAD with the 'corrected' dose/time offset
assignment, unfiltered and after the connectivity-score Monte Carlo p-value
filter, over the per-instance CS runs of tcga_wsRGES_cs_runs_used.tsv.

Design
------
- Instances: the CS runs listed in OUT_DIR/tcga_wsRGES_cs_runs_used.tsv
  (cancer x configuration x LINCS cell line), loaded with loader.load_cs_run
  (selected_cs_run_id); instance metadata from conf.LINCS_METADATA_PATH;
  cell-line weights: the 64 lines of OUT_DIR/tcga_cell_line_weights.tsv (TCGA
  tumour correlation of the DepMap, HPA and GEO profiles) in every
  configuration, offset assignment and filter variant; the 58 DepMap lines of
  OUT_DIR/tcga_cell_line_weights_ccle58.tsv (the Chen 2017 line set) only in
  the Chen 2017 reproduction ('binchen_r', unfiltered); plus
  loader.cor_weights_for_disease (Chen 2017 correlation).
- Collapses of wsRGES_IC50_correlations.collapse (drug = pert_iname):
  'wsRGES', 'wsRGES_chen2017w', 'sRGES_all', 'best_all' over all lines, and
  'best' / 'srges' over the validation line of the cancer (MCF7, HEPG2, HT29).
  'best' and 'best_all' use no offsets: identical in every assignment.
  Only the dose/time offset assignment (sRGES.compute_sRGES, argument
  srges_diff_mode) differs from that function: here a parameter.
- Unfiltered: assignment 'corrected' (each instance receives the offset of its
  own dose bin and time). The same code with 'binchen_r' and the DepMap
  weights only (OUT_DIR/tcga_cell_line_weights_ccle58.tsv, the weights of
  that table) reproduces OUT_DIR/tcga_wsRGES_per_drug.tsv; the largest absolute score difference and
  the Spearman correlation between the two assignments (both with the
  64-line weights) are written to
  OUT_DIR/tcga_wsRGES_corrected_offsets_summary.tsv. The 'DEG' rows of that
  reproduction are the Chen 2017 reproduction table
  OUT_DIR/tcga_wsRGES_per_drug_DEG__binchen_r.tsv.
- Filter (src/case_study/cs_pvalue_filter.flag): p' = (k + 1) / 1001 with
  k = cs_p_value x 1000, Benjamini-Hochberg adjustment within one CS run (one
  cancer x configuration x cell line), instances with p_adj <= 0.05 are kept.
  The collapses are recomputed from the kept instances with 'corrected':
  'bh05' rebuilds the offset table from the kept instances, 'bh05f' (frozen
  offsets) scores the kept instances with the offset table of the unfiltered
  instances (all lines for the all-line collapses, the validation line for
  'srges'). A drug without a kept instance leaves the table.
- Sign convention: more negative = stronger reversal of the disease signature.

Inputs
------
- OUT_DIR/tcga_wsRGES_cs_runs_used.tsv, OUT_DIR/tcga_cell_line_weights.tsv,
  OUT_DIR/tcga_cell_line_weights_ccle58.tsv, OUT_DIR/tcga_wsRGES_per_drug.tsv
- conf.LINCS_METADATA_PATH, the CS runs of loader.load_cs_run

Outputs (OUT_DIR = conf.BINCHEN_VALIDATION_OUT_DIR)
-------
- OUT_DIR/tcga_wsRGES_per_drug_DEG__binchen_r.tsv
    configuration 'DEG', assignment 'binchen_r', 58 lines.
- OUT_DIR/tcga_wsRGES_per_drug_corrected.tsv
    columns of tcga_wsRGES_per_drug.tsv.
- OUT_DIR/tcga_wsRGES_per_drug_corrected_pvalue_filtered.tsv
- OUT_DIR/tcga_wsRGES_per_drug_corrected_pvalue_filtered_frozen_offsets.tsv
    same columns.
- OUT_DIR/tcga_wsRGES_corrected_offsets_summary.tsv
    one row per (cancer, config, collapse): recomputation check of the
    'binchen_r' assignment (max |score difference|, drugs compared), Spearman
    correlation of the 'binchen_r' and 'corrected' scores on the shared drugs,
    and for the p-value filter the instances before / after and the drugs
    before / after.

Usage
-----
    cd src && PYTHONPATH=. python validations/chembl/wsRGES_tcga_corrected_offsets.py
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import conf
import sRGES
from loader import load_cs_run, read_lincs_metadata, cor_weights_for_disease
from sRGES import compute_sRGES
from case_study import cs_pvalue_filter as pf
from case_study import per_drug_scores as pds
from validations.chembl.final_step_correlations_IC50_CS import collapse_profiles_to_drug
from validations.chembl.wsRGES_IC50_correlations import CANCERS, DRUG, LINE_COLLAPSES

OUT_DIR = conf.BINCHEN_VALIDATION_OUT_DIR
COLS = ['cancer', 'config', 'collapse', DRUG, 'score', 'n_instances', 'n_cell_lines']
CORR_MODE = 'corrected'


def collapse(inst, md, line, w_tcga, w_chen, mode, frozen_all=None, frozen_line=None):
    """Collapses of wsRGES_IC50_correlations.collapse with srges_diff_mode = `mode`;
    frozen_all / frozen_line: offset tables replacing the rebuilt ones (all-line / validation-line collapses)."""
    x = inst.assign(LINCS_id=inst['LINCS_id'].astype(str)).merge(md, left_on='LINCS_id', right_on='id')
    for c in ('pert_dose', 'pert_time', 'connectivity_score'):
        x[c] = pd.to_numeric(x[c], errors='coerce')
    kept = x[(x['pert_dose'] > 0) & x['pert_time'].isin([6, 24])]
    out = []

    def run(fn, frozen):
        if frozen is None:
            return fn()
        with pf.frozen_offsets(frozen):
            return fn()

    for name, w in [('wsRGES', w_tcga), ('wsRGES_chen2017w', w_chen), ('sRGES_all', None)]:
        s = run(lambda: compute_sRGES(x, score_col='connectivity_score', drug_col=DRUG, cell_col='cell_id',
                                      cell_line_weights=w, srges_diff_mode=mode), frozen_all)
        k = kept[kept['cell_id'].isin(w)] if w is not None else kept
        out.append(pd.DataFrame({'collapse': name, DRUG: s[DRUG], 'score': s['sRGES'], 'n_instances': s['n'],
                                 'n_cell_lines': s[DRUG].map(k.groupby(DRUG)['cell_id'].nunique())}))
    g = kept.groupby(DRUG)
    out.append(pd.DataFrame({'collapse': 'best_all', DRUG: g['connectivity_score'].min().index,
                             'score': g['connectivity_score'].min().values, 'n_instances': g.size().values,
                             'n_cell_lines': g['cell_id'].nunique().values}))
    xl = x[x['cell_id'] == line][['LINCS_id', DRUG, 'connectivity_score']]
    for how in LINE_COLLAPSES:
        kw = {'srges_kwargs': {'srges_diff_mode': mode}} if how == 'srges' else {}
        s = run(lambda: collapse_profiles_to_drug(xl, score_col='connectivity_score', drug_col=DRUG, how=how, **kw),
                frozen_line if how == 'srges' else None)
        out.append(pd.DataFrame({'collapse': how, DRUG: s[DRUG], 'score': s['connectivity_score'],
                                 'n_instances': s[DRUG].map(xl.groupby(DRUG).size()), 'n_cell_lines': 1}))
    return pd.concat(out, ignore_index=True)


def load_instances(cancer, label, used):
    parts = []
    for _, r in used[(used['cancer'] == cancer) & (used['config'] == label)].iterrows():
        df, _, run_id = load_cs_run(disease=cancer, landmark_disease=False, selected_cs_run_id=r['cs_run_id'],
                                    verbose=False)
        assert run_id == r['cs_run_id'] and len(df) == r['n_rows']
        df = df[['LINCS_id', 'connectivity_score', 'cs_p_value']].reset_index(drop=True)
        df = pd.concat([df, pf.flag(df['cs_p_value'])], axis=1).assign(cell_line=r['cell_line'])
        parts.append(df)
    return pd.concat(parts, ignore_index=True)


def offsets_of(inst, md):
    return pds.offsets(inst[['LINCS_id', 'connectivity_score']], md)


if __name__ == '__main__':
    assert conf.srges_diff_mode == 'corrected'
    md = read_lincs_metadata(dtype={'id': str})[['id', 'pert_id', DRUG, 'cell_id', 'pert_dose', 'pert_time']]
    wt = pd.read_csv(OUT_DIR / 'tcga_cell_line_weights.tsv', sep='\t')
    wt58 = pd.read_csv(OUT_DIR / 'tcga_cell_line_weights_ccle58.tsv', sep='\t')
    used = pd.read_csv(OUT_DIR / 'tcga_wsRGES_cs_runs_used.tsv', sep='\t')
    old = pd.read_csv(OUT_DIR / 'tcga_wsRGES_per_drug.tsv', sep='\t')
    tabs = {m: [] for m in ('binchen_r', 'binchen_r58', 'corrected', 'bh05', 'bh05f')}
    summ = []
    for cancer, line in CANCERS.items():
        w_tcga = dict(zip(wt.loc[wt.cancer == cancer, 'lincs_cell_id'], wt.loc[wt.cancer == cancer, 'cor']))
        w_58 = dict(zip(wt58.loc[wt58.cancer == cancer, 'lincs_cell_id'], wt58.loc[wt58.cancer == cancer, 'cor']))
        w_chen = cor_weights_for_disease(cancer)
        for label in used.loc[used['cancer'] == cancer, 'config'].unique():
            inst = load_instances(cancer, label, used)
            print(f'{cancer} {label}: {inst["cell_line"].nunique()} lines, {len(inst)} instances', flush=True)
            keep = inst['keep_bh05']
            f_all = offsets_of(inst, md)
            f_line = offsets_of(inst[inst['LINCS_id'].astype(str).isin(
                md.loc[md['cell_id'] == line, 'id'])], md)
            res = {
                'binchen_r58': collapse(inst, md, line, w_58, w_chen, 'binchen_r'),
                'binchen_r': collapse(inst, md, line, w_tcga, w_chen, 'binchen_r'),
                'corrected': collapse(inst, md, line, w_tcga, w_chen, CORR_MODE),
                'bh05': collapse(inst[keep], md, line, w_tcga, w_chen, CORR_MODE),
                'bh05f': collapse(inst[keep], md, line, w_tcga, w_chen, CORR_MODE, f_all, f_line)}
            for m, t in res.items():
                tabs[m].append(t.assign(cancer=cancer, config=label)[COLS])
            o = old[(old['cancer'] == cancer) & (old['config'] == label)]
            for coll, a in res['corrected'].groupby('collapse'):
                b = res['binchen_r'][res['binchen_r']['collapse'] == coll].set_index(DRUG)['score']
                b58 = res['binchen_r58'][res['binchen_r58']['collapse'] == coll].set_index(DRUG)['score']
                oo = o[o['collapse'] == coll].set_index(DRUG)['score']
                a = a.set_index(DRUG)['score']
                sh = a.index.intersection(b.index)
                j = b58.to_frame('new').join(oo.rename('old'), how='outer')
                fl = res['bh05'][res['bh05']['collapse'] == coll]
                summ.append({'cancer': cancer, 'config': label, 'collapse': coll,
                             'n_drugs_published': len(oo), 'n_drugs_recomputed': len(b58),
                             'max_abs_diff_binchen_r_vs_stored': float((j['new'] - j['old']).abs().max()),
                             'n_drugs_only_one_table': int(j['new'].isna().sum() + j['old'].isna().sum()),
                             'spearman_binchen_r_vs_corrected': spearmanr(a[sh], b[sh]).statistic,
                             'n_instances_before': len(inst), 'n_instances_kept': int(keep.sum()),
                             'n_drugs_before': len(a), 'n_drugs_after_filter': len(fl)})
    for m, fn in [('corrected', 'tcga_wsRGES_per_drug_corrected.tsv'),
                  ('bh05', 'tcga_wsRGES_per_drug_corrected_pvalue_filtered.tsv'),
                  ('bh05f', 'tcga_wsRGES_per_drug_corrected_pvalue_filtered_frozen_offsets.tsv')]:
        pd.concat(tabs[m], ignore_index=True).to_csv(OUT_DIR / fn, sep='\t', index=False)
    deg = pd.concat(tabs['binchen_r58'], ignore_index=True)
    deg[deg['config'] == 'DEG'].to_csv(OUT_DIR / 'tcga_wsRGES_per_drug_DEG__binchen_r.tsv', sep='\t', index=False)
    pd.DataFrame(summ).to_csv(OUT_DIR / 'tcga_wsRGES_corrected_offsets_summary.tsv', sep='\t', index=False)
    s = pd.DataFrame(summ)
    print(s.groupby('collapse')[['max_abs_diff_binchen_r_vs_stored', 'n_drugs_only_one_table']].max())
