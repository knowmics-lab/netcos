#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-drug summarized scores of BRCA, LIHC and COAD over the per-instance
connectivity scores of all LINCS cell lines, and their correlation with the
ChEMBL IC50 of the Chen 2017 validation lines (MCF7, HepG2, HT29).

Design
------
- Per-instance CS: for every LINCS cell line, the most recent cs_batch run of
  the cancer (disease_run_id = cancer, Chen 2017 TCGA signature) on drug input
  '<cell_id>_LM', CS_METHOD = 'bin_chen', resolved with loader.load_cs_run
  (drug_run_id fixed, so a run of one line never stands in for another). DEG
  is scored only when every line has a run; a MITHrIL configuration is scored
  on the lines that have a run (a line whose MITHrIL inputs are incomplete
  never gets one), and the lines left out are printed; the number of runs
  used is in the appended registry rows.
- Configurations: CONFIGS (DEG and the three MITHrIL configurations). A
  MITHrIL configuration without any per-instance run is reported and skipped;
  its collapses are named '<config>_<collapse>' (e.g. 'mith_wsRGES').
- Per-drug collapses (drug = pert_iname, as in final_step_correlations_IC50_CS):
  * 'wsRGES': sRGES.compute_sRGES over the instances of all lines, each line
    weighted by its median Spearman correlation with the cancer's TCGA primary
    tumours (tcga_ccle_cell_line_weights.py, cor / max(cor)); lines without
    DepMap expression are left out.
  * 'wsRGES_chen2017w': as 'wsRGES' with the Chen 2017 per-line correlations
    (cell_line_<T>_tacle.csv, loader.cor_weights_for_disease).
  * 'sRGES_all': compute_sRGES over all lines, every weight = 1.
  * 'best_all': the most negative CS among the drug's instances with
    pert_dose > 0 and pert_time in {6, 24} h, all lines.
  * 'best' and 'srges' of the validation line only, computed exactly as the
    'best' / 'srges' LINCS collapses of final_step_correlations_IC50_CS
    (collapse_profiles_to_drug on that line's instances).
  compute_sRGES uses the Chen 2017 dose/time correction with the index mapping
  of the published R code (srges_diff_mode 'binchen_r').
- IC50 side, the same grid as the existing validation: SD8 IC50 rows only
  (ic50_only), collapsed per pert_iname by 'best' or 'median', restricted to
  the validation line ('matched') or all lines of the SD8 sheet ('all').
  Two drug sets per all-line collapse: 'all_scored' (every drug with a
  score) and 'line_drugs' (drugs with at least one instance on the
  validation line, the drug set of 'best' / 'srges'); with the 'best' and
  'srges' collapses the two sets coincide and only 'all_scored' is reported.
  Merge with _merge_cs_ic50_inchi_first; Spearman of score vs IC50 (linear and
  log10, identical by construction), Pearson of score vs log10 IC50,
  precision/recall from classify_ic50_vs_cs with the thresholds of the
  existing grid (CS_TH, IC50_TH; IC50 in SD8 units).
- Reference values: the rows of conf.chembl_val_log_filename computed on
  the validation line alone (drug_run_id '<line>_LM', CS_METHOD 'bin_chen',
  CS runs not before CS_CUTOFF): the DEG rows (landmark_disease False,
  cs_on_LM 0, LINCS collapses 'best' / 'srges'), the configuration with the
  highest Spearman per IC50 recipe, and Chen 2017 SD5 sRGES vs IC50
  ('spearman_r_SD5').
- Paired bootstrap of Spearman differences: for every cancer, IC50 recipe
  and pair of collapses (line_drugs set; 'best' / 'srges' on their own drugs,
  the same set), the drugs are resampled with replacement N_BOOT times
  (seed BOOT_SEED) and the difference of the two Spearman correlations with
  IC50 is recomputed on each resample; percentile 95% interval. Intervals are
  not adjusted for the number of comparisons (6 collapses x 4 IC50 recipes x
  3 cancers).
- Agreement of each per-drug score with the Chen 2017 published sRGES of the
  cancer (lincs_cancer_sRGES.csv), Spearman over shared pert_iname.
- Sign convention: more negative = stronger reversal; a positive Spearman
  with IC50 means stronger reversal goes with lower IC50.

Inputs
------
- conf.cs_log_filename and the CS files it points to
- conf.LINCS_METADATA_PATH, conf.ic50_file (SD8), conf.BC_DISEASE_DATA
- OUT_DIR/tcga_cell_line_weights.tsv
- conf.chembl_val_log_filename (existing grid rows)

Outputs (OUT_DIR = conf.BINCHEN_VALIDATION_OUT_DIR, IMG_OUT = conf.IMG_DIR / 'binchen2017_validation')
-------
- OUT_DIR/tcga_wsRGES_per_drug.tsv
    one row per (cancer, config, collapse, pert_iname): score, n_instances,
    n_cell_lines.
- OUT_DIR/tcga_wsRGES_cs_runs_used.tsv
    cs_run_id used for every (cancer, config, cell line).
- OUT_DIR/tcga_wsRGES_vs_IC50.tsv
    one row per (cancer, config, collapse, drug_set, IC50 recipe): n_drugs, spearman_r,
    spearman_pval, pearson_r_log, pearson_pval_log, tp, fp, fn, tn,
    precision, recall, spearman_vs_chen2017_sRGES, n_vs_chen2017; plus the
    reference rows (source 'grid_DEG', 'grid_best_config', 'chen2017_SD5').
- OUT_DIR/tcga_wsRGES_bootstrap_spearman_diff.tsv
    one row per (cancer, IC50 recipe, collapse_a, collapse_b), DEG,
    line_drugs: n_drugs, rho_a, rho_b, diff = rho_a - rho_b, ci_low,
    ci_high, frac_boot_diff_le0, n_boot.
- conf.chembl_val_log_filename (only when LOG_APPEND): one appended row per
  (cancer, config, all-line collapse, drug set, IC50 recipe), schema of the
  existing rows;
  LINCS_drug_collapse_method = collapse name ('<collapse>_lineDrugs' for the
  'line_drugs' set), drug_run_id = 'all_LM',
  cs_run_id = first and last run ids and the number of runs, cs_file =
  OUT_DIR/tcga_wsRGES_cs_runs_used.tsv.
- IMG_OUT/tcga_wsRGES_vs_IC50.png
    top: wsRGES vs log10 IC50 (matched line, median IC50) per cancer, drugs
    with instances on the validation line, with the bootstrap intervals of
    wsRGES - srges and wsRGES - best; bottom: Spearman with IC50 of every
    collapse for the four IC50 recipes, same drugs.

Usage
-----
    cd src && PYTHONPATH=. python validations/chembl/wsRGES_IC50_correlations.py
"""

import warnings
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, rankdata

warnings.filterwarnings('ignore', message='Unknown extension is not supported', category=UserWarning)

import conf
import lfs_style
from loader import load_cs_run, load_IC50, cor_weights_for_disease, read_lincs_metadata
from logger import append_run_metadata
from sRGES import compute_sRGES
from validations.chembl.final_step_correlations_IC50_CS import (
    collapse_profiles_to_drug, _merge_cs_ic50_inchi_first, classify_ic50_vs_cs)

OUT_DIR = conf.BINCHEN_VALIDATION_OUT_DIR
IMG_OUT = conf.IMG_DIR / 'binchen2017_validation'
CANCERS = {'BRCA': 'MCF7', 'LIHC': 'HEPG2', 'COAD': 'HT29'}
CONFIGS = {  # label: (mith, cs_on_LM, CS_ON_PATHWAYS)
    'DEG': (0, 0, 0),
    'mith': (1, 0, 0),
    'mith_LM': (1, 1, 0),
    'mith_pw': (1, 0, 1),
}
RUN_SUFFIX = '_LM'
DRUG = 'pert_iname'
ALL_LINE_COLLAPSES = ['wsRGES', 'wsRGES_chen2017w', 'sRGES_all', 'best_all']
LINE_COLLAPSES = ['best', 'srges']
IC50_COLLAPSES = ['best', 'median']
CS_TH, IC50_TH = -1.5, 10.0  # thresholds of the existing grid
CS_CUTOFF = pd.Timestamp('2026-05-11')  # bin_chen runs before this date are not used
LOG_APPEND = False  # True appends the all-line rows to conf.chembl_val_log_filename
N_BOOT, BOOT_SEED = 5000, 0
BOOT_ORDER = ['wsRGES', 'wsRGES_chen2017w', 'sRGES_all', 'best_all', 'srges', 'best']
DRUG_SETS = ['all_scored', 'line_drugs']
REQUIRE_ALL_LINES = {'DEG'}  # MITHrIL configurations are scored on the lines that have runs


def load_instances(cancer, label, cells):
    m, lm, pw = CONFIGS[label]
    parts, used, missing = [], [], []
    for c in cells:
        try:
            df, path, run_id = load_cs_run(disease=cancer, landmark_disease=False, drug_run_id=c + RUN_SUFFIX,
                                           mith=m, cs_on_LM=lm, CS_ON_PATHWAYS=pw, CS_METHOD='bin_chen',
                                           verbose=False)
        except (FileNotFoundError, ValueError):
            missing.append(c)
            continue
        if pd.to_datetime(run_id[:16], format='%d_%m_%Y_%H_%M') < CS_CUTOFF:
            missing.append(c)
            continue
        parts.append(df[['LINCS_id', 'connectivity_score']])
        used.append({'cancer': cancer, 'config': label, 'cell_line': c, 'cs_run_id': run_id, 'n_rows': len(df)})
    inst = pd.concat(parts, ignore_index=True) if parts else None
    return inst, used, missing


def collapse(inst, md, line, w_tcga, w_chen):
    x = inst.assign(LINCS_id=inst['LINCS_id'].astype(str)).merge(md, left_on='LINCS_id', right_on='id')
    for c in ('pert_dose', 'pert_time', 'connectivity_score'):
        x[c] = pd.to_numeric(x[c], errors='coerce')
    kept = x[(x['pert_dose'] > 0) & x['pert_time'].isin([6, 24])]
    out = []
    for name, w in [('wsRGES', w_tcga), ('wsRGES_chen2017w', w_chen), ('sRGES_all', None)]:
        s = compute_sRGES(x, score_col='connectivity_score', drug_col=DRUG, cell_col='cell_id',
                          cell_line_weights=w, srges_diff_mode='binchen_r')
        k = kept[kept['cell_id'].isin(w)] if w is not None else kept
        out.append(pd.DataFrame({'collapse': name, DRUG: s[DRUG], 'score': s['sRGES'], 'n_instances': s['n'],
                                 'n_cell_lines': s[DRUG].map(k.groupby(DRUG)['cell_id'].nunique())}))
    g = kept.groupby(DRUG)
    out.append(pd.DataFrame({'collapse': 'best_all', DRUG: g['connectivity_score'].min().index,
                             'score': g['connectivity_score'].min().values, 'n_instances': g.size().values,
                             'n_cell_lines': g['cell_id'].nunique().values}))
    xl = x[x['cell_id'] == line][['LINCS_id', DRUG, 'connectivity_score']]
    for how in LINE_COLLAPSES:
        s = collapse_profiles_to_drug(xl, score_col='connectivity_score', drug_col=DRUG, how=how)
        out.append(pd.DataFrame({'collapse': how, DRUG: s[DRUG], 'score': s['connectivity_score'],
                                 'n_instances': s[DRUG].map(xl.groupby(DRUG).size()), 'n_cell_lines': 1}))
    return pd.concat(out, ignore_index=True)


def load_ic50_grid(cancer, line):
    grid = {}
    for cl in (line, None):
        raw, log = load_IC50(conf.ic50_file, cancer, cell_line=cl, IC50_ONLY=True, IC_50_binchen_SD5=True)
        for how in IC50_COLLAPSES:
            keep = ['standard_inchi'] if 'standard_inchi' in raw.columns else []
            ic = collapse_profiles_to_drug(raw, score_col='standard_value', drug_col=DRUG, how=how, keep_cols=keep)
            grid[('matched' if cl else 'all', how)] = (ic, len(raw), log)
    return grid


def correlate(scores, ic50):
    dr = scores.rename(columns={'score': 'connectivity_score'})[[DRUG, 'connectivity_score']]
    m = _merge_cs_ic50_inchi_first(dr, ic50, cs_drug_colname=DRUG, ic50_drug_colname=DRUG)
    m = m.dropna(subset=['connectivity_score', 'standard_value'])
    rho, p = spearmanr(m['connectivity_score'], m['standard_value'])
    pr, pp = pearsonr(m['connectivity_score'], np.log10(m['standard_value']))
    _, met = classify_ic50_vs_cs(m, score_col='connectivity_score', ic50_col='standard_value',
                                 cs_threshold=CS_TH, ic50_threshold=IC50_TH)
    return m, {'n_drugs': m[DRUG].nunique(), 'spearman_r': rho, 'spearman_pval': p,
               'pearson_r_log': pr, 'pearson_pval_log': pp,
               **{k: met[k] for k in ('tp', 'fp', 'fn', 'tn', 'precision', 'recall')}}


def grid_reference():
    r = pd.read_csv(conf.chembl_val_log_filename, sep='\t')
    r = r[r['drug_run_id'].astype(str) == r['cell_line'].astype(str) + RUN_SUFFIX].copy()
    r = r[r['cs_run_id'].astype(str).str.fullmatch(r'\d\d_\d\d_\d{4}_\d\d_\d\d_.*_bin_chen_connectivity_score')]
    r = r[pd.to_datetime(r['cs_run_id'].str[:16], format='%d_%m_%Y_%H_%M') >= CS_CUTOFF]
    r['ic50_recipe_cells'] = np.where(r['cell_line_IC50'].isna(), 'all', 'matched')
    r['config'] = r['cs_run_id'].str.extract(r'\d\d_\d\d_\d{4}_\d\d_\d\d_(.*)_bin_chen_connectivity_score')[0]
    r['cancer'] = r['cell_line'].map({v: k for k, v in CANCERS.items()})
    r = r.drop_duplicates(['cs_run_id', 'LINCS_drug_collapse_method', 'IC50_drug_collapse_method',
                           'ic50_recipe_cells', 'landmark_disease', 'cs_threshold'], keep='last')
    sd5_all = r.dropna(subset=['spearman_r_SD5'])
    r = r[r['cs_threshold'] == CS_TH]
    rows = []
    deg = r[(r['config'] == 'DEG') & (~r['landmark_disease'].astype(str).str.lower().isin(['true', '1']))]
    for _, x in deg.iterrows():
        rows.append({'source': 'grid_DEG', 'cancer': x.cancer, 'config': 'DEG',
                     'collapse': x.LINCS_drug_collapse_method,
                     'ic50_cells': x.ic50_recipe_cells, 'ic50_collapse': x.IC50_drug_collapse_method,
                     'n_drugs': x.n_unique_drugs, 'spearman_r': x.spearman_r, 'cs_run_id': x.cs_run_id})
    for (cancer, cells, ich), g in r.groupby(['cancer', 'ic50_recipe_cells', 'IC50_drug_collapse_method']):
        b = g.loc[g['spearman_r'].idxmax()]
        lm = '_LMdisease' if str(b.landmark_disease).lower() in ('true', '1') else ''
        rows.append({'source': 'grid_best_config', 'cancer': cancer, 'config': f'{b.config}{lm}',
                     'collapse': b.LINCS_drug_collapse_method, 'ic50_cells': cells, 'ic50_collapse': ich,
                     'n_drugs': b.n_unique_drugs, 'spearman_r': b.spearman_r, 'cs_run_id': b.cs_run_id})
    for cancer, g in sd5_all.groupby('cancer'):
        rows.append({'source': 'chen2017_SD5', 'cancer': cancer, 'config': 'Chen2017', 'collapse': 'sRGES',
                     'ic50_cells': 'matched', 'ic50_collapse': 'SD5', 'n_drugs': g['n_overlap_BC_SD5'].iloc[-1],
                     'spearman_r': g['spearman_r_SD5'].iloc[-1]})
    return pd.DataFrame(rows).drop_duplicates(['source', 'cancer', 'collapse', 'ic50_cells', 'ic50_collapse'])


def spearman_cols(scores, y):
    """Spearman of every column of scores (n x k) with y (n), average ranks."""
    rs = np.apply_along_axis(rankdata, 0, scores)
    ry = rankdata(y)
    rs = rs - rs.mean(axis=0)
    ry = ry - ry.mean()
    return (rs * ry[:, None]).sum(axis=0) / np.sqrt((rs ** 2).sum(axis=0) * (ry ** 2).sum())


def bootstrap_diffs(boot_in):
    rows = []
    rng = np.random.default_rng(BOOT_SEED)
    for (cancer, cells, ich), frames in boot_in.items():
        cols = [c for c in BOOT_ORDER if c in frames]
        d = pd.concat([frames[c]['connectivity_score'].rename(c) for c in cols], axis=1, join='inner')
        y = frames[cols[0]]['standard_value'].reindex(d.index).to_numpy(float)
        X = d.to_numpy(float)
        rho = spearman_cols(X, y)
        bs = np.array([spearman_cols(X[i], y[i]) for i in rng.integers(0, len(y), (N_BOOT, len(y)))])
        for a in range(len(cols)):
            for b in range(a + 1, len(cols)):
                diff = bs[:, a] - bs[:, b]
                lo, hi = np.nanpercentile(diff, [2.5, 97.5])
                rows.append({'cancer': cancer, 'ic50_cells': cells, 'ic50_collapse': ich, 'n_drugs': len(y),
                             'collapse_a': cols[a], 'collapse_b': cols[b], 'rho_a': rho[a], 'rho_b': rho[b],
                             'diff': rho[a] - rho[b], 'ci_low': lo, 'ci_high': hi,
                             'frac_boot_diff_le0': float(np.mean(diff <= 0)), 'n_boot': N_BOOT})
    return pd.DataFrame(rows)


def log_row(cancer, line, label, coll, cells, ich, res, used, n_ic50, ic_log):
    ids = sorted(u['cs_run_id'] for u in used)
    row = {
        'IC50_validation_run_id': datetime.now().strftime('%d_%m_%Y_%H_%M_%S'),
        'datetime': datetime.now().isoformat(),
        'cs_file': str(OUT_DIR / 'tcga_wsRGES_cs_runs_used.tsv'),
        'ic50_file': str(conf.ic50_file), 'lincs_metadata_file': str(conf.LINCS_METADATA_PATH),
        'disease_run_id': cancer, 'drug_run_id': 'all' + RUN_SUFFIX,
        'cs_drug_colname': DRUG, 'ic50_drug_colname': DRUG,
        'LINCS_drug_collapse_method': coll if label == 'DEG' else f'{label}_{coll}',
        'IC50_drug_collapse_method': ich,
        'n_ic50_rows_before_drug_collapse': n_ic50,
        'n_merged_rows': res['n_drugs'], 'n_unique_drugs': res['n_drugs'],
        'spearman_r': round(res['spearman_r'], 2), 'spearman_pval': round(res['spearman_pval'], 2),
        'spearman_r_log': round(res['spearman_r'], 2), 'spearman_pval_log': round(res['spearman_pval'], 2),
        **{k: res[k] for k in ('tp', 'fp', 'fn', 'tn')},
        'precision': res['precision'], 'recall': res['recall'],
        'cs_run_id': f'{ids[0]} .. {ids[-1]} ({len(ids)} runs)',
        'cell_line': line, 'landmark_disease': False, 'landmark_drug': True, 'CS_METHOD': 'bin_chen',
        'cs_mith': CONFIGS[label][0], 'cs_on_LM': CONFIGS[label][1], 'CS_ON_PATHWAYS': bool(CONFIGS[label][2]),
        'ic50_only': True, 'cell_line_IC50': line if cells == 'matched' else None,
        'cs_threshold': CS_TH, 'ic50_threshold': IC50_TH,
    }
    row.update(ic_log)
    append_run_metadata(conf.chembl_val_log_filename, row)


def figure(res, per_drug_merged, boot):
    lfs_style.apply()
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    order = ['best', 'srges', 'best_all', 'sRGES_all', 'wsRGES', 'wsRGES_chen2017w']
    recipes = [('matched', 'median'), ('matched', 'best'), ('all', 'median'), ('all', 'best')]
    for j, (cancer, line) in enumerate(CANCERS.items()):
        a = ax[0, j]
        m = per_drug_merged[cancer]
        a.scatter(m['connectivity_score'], np.log10(m['standard_value']), s=14, facecolors='none',
                  edgecolors=lfs_style.BEFORE, linewidths=.7)
        r = res[(res.cancer == cancer) & (res.config == 'DEG') & (res.collapse == 'wsRGES')
                & (res.drug_set == 'line_drugs') & (res.ic50_cells == 'matched')
                & (res.ic50_collapse == 'median')].iloc[0]
        a.set_xlabel('wsRGES (DEG)')
        a.set_ylabel(f'log10 IC50, {line} (median)')
        a.set_title(f'{cancer}, drugs with instances on {line}: n={int(r.n_drugs)}, '
                    f'Spearman {r.spearman_r:.2f}, Pearson {r.pearson_r_log:.2f}', fontsize=10)
        b = boot[(boot.cancer == cancer) & (boot.ic50_cells == 'matched') & (boot.ic50_collapse == 'median')
                 & (boot.collapse_a == 'wsRGES') & boot.collapse_b.isin(LINE_COLLAPSES)].set_index('collapse_b')
        txt = '\n'.join(f'wsRGES - {c}: {b.loc[c, "diff"]:+.2f} [{b.loc[c, "ci_low"]:+.2f}, {b.loc[c, "ci_high"]:+.2f}]'
                        for c in ('srges', 'best'))
        a.text(0.02, 0.98, txt + f'\n(paired bootstrap, {N_BOOT} resamples, 95%)', transform=a.transAxes,
               va='top', ha='left', fontsize=8, bbox=dict(boxstyle='round', facecolor='white', alpha=.85))
        a = ax[1, j]
        sub = res[(res.cancer == cancer) & (res.config == 'DEG')
                  & ((res.drug_set == 'line_drugs') | res.collapse.isin(LINE_COLLAPSES))]
        width = .2
        for k, (cells, ich) in enumerate(recipes):
            s = sub[(sub.ic50_cells == cells) & (sub.ic50_collapse == ich)].set_index('collapse').reindex(order)
            a.bar(np.arange(len(order)) + (k - 1.5) * width, s['spearman_r'], width,
                  label=f'IC50 {cells} line(s), {ich}', color=['#0000FF', '#7f7fff', '#FF0000', '#ff9f9f'][k])
        a.axhline(0, color='k', lw=.6)
        a.set_xticks(np.arange(len(order)))
        a.set_xticklabels(order, rotation=30, fontsize=8)
        a.set_ylabel('Spearman (score vs IC50)')
        a.set_title(f'{cancer}: drugs with instances on {line}', fontsize=10)
    h, l = ax[1, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=4, fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    IMG_OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMG_OUT / 'tcga_wsRGES_vs_IC50.png', dpi=150)


if __name__ == '__main__':
    md = read_lincs_metadata(dtype={'id': str})
    md = md[['id', 'pert_id', DRUG, 'cell_id', 'pert_dose', 'pert_time']]
    cells = sorted(md['cell_id'].unique())
    wt = pd.read_csv(OUT_DIR / 'tcga_cell_line_weights.tsv', sep='\t')
    print(f'{len(cells)} LINCS cell lines')

    per_drug, used_all, res_rows, merged_fig, boot_in = [], [], [], {}, {}
    for cancer, line in CANCERS.items():
        w_tcga = dict(zip(wt.loc[wt.cancer == cancer, 'lincs_cell_id'], wt.loc[wt.cancer == cancer, 'cor']))
        w_chen = cor_weights_for_disease(cancer)
        print(f'\n{cancer}: weights TCGA/DepMap {len(w_tcga)} lines, Chen 2017 {len(w_chen)} lines')
        ic_grid = load_ic50_grid(cancer, line)
        chen = pd.read_csv(conf.BC_DISEASE_DATA / cancer / 'lincs_cancer_sRGES.csv')[[DRUG, 'sRGES']]
        chen = chen.groupby(DRUG)['sRGES'].mean()
        for label in CONFIGS:
            inst, used, missing = load_instances(cancer, label, cells)
            if inst is None or (missing and label in REQUIRE_ALL_LINES):
                print(f'{cancer} {label}: per-instance CS missing for {len(missing)}/{len(cells)} lines; '
                      f'{label}_wsRGES / {label}_sRGES_all not computed')
                continue
            if missing:
                print(f'{cancer} {label}: scored on {len(used)} lines; no per-instance CS for '
                      f'{len(missing)} lines: {", ".join(missing)}')
            used_all += used
            print(f'{cancer} {label}: {len(used)} lines, {len(inst)} instances, runs '
                  f'{min(u["cs_run_id"] for u in used)} .. {max(u["cs_run_id"] for u in used)}')
            sc = collapse(inst, md, line, w_tcga, w_chen)
            sc.insert(0, 'config', label)
            sc.insert(0, 'cancer', cancer)
            per_drug.append(sc)
            line_drugs = set(sc.loc[sc.collapse == 'best', DRUG])
            for coll, s_all in sc.groupby('collapse'):
                v = s_all.set_index(DRUG)['score'].to_frame().join(chen, how='inner')
                rho_chen = spearmanr(v['score'], v['sRGES']).statistic
                for drug_set in DRUG_SETS:
                    if drug_set == 'line_drugs' and coll in LINE_COLLAPSES:
                        continue
                    s = s_all[s_all[DRUG].isin(line_drugs)] if drug_set == 'line_drugs' else s_all
                    for (cells_, ich), (ic, n_ic, ic_log) in ic_grid.items():
                        m, r = correlate(s, ic)
                        res_rows.append({'source': 'this_run', 'cancer': cancer, 'config': label,
                                         'collapse': coll, 'drug_set': drug_set,
                                         'ic50_cells': cells_, 'ic50_collapse': ich, **r,
                                         'spearman_vs_chen2017_sRGES': rho_chen, 'n_vs_chen2017': len(v)})
                        on_line = drug_set == 'line_drugs' or coll in LINE_COLLAPSES
                        if label == 'DEG' and on_line:
                            boot_in.setdefault((cancer, cells_, ich), {})[coll] = \
                                m.drop_duplicates(DRUG).set_index(DRUG)[['connectivity_score', 'standard_value']]
                        if (label == 'DEG' and coll == 'wsRGES' and drug_set == 'line_drugs'
                                and cells_ == 'matched' and ich == 'median'):
                            merged_fig[cancer] = m
                        if LOG_APPEND and coll in ALL_LINE_COLLAPSES:
                            name = coll if drug_set == 'all_scored' else f'{coll}_lineDrugs'
                            log_row(cancer, line, label, name, cells_, ich, r, used, n_ic, ic_log)

    res = pd.DataFrame(res_rows)
    ref = grid_reference()
    out = pd.concat([res, ref], ignore_index=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.concat(per_drug, ignore_index=True).to_csv(OUT_DIR / 'tcga_wsRGES_per_drug.tsv', sep='\t', index=False)
    pd.DataFrame(used_all).to_csv(OUT_DIR / 'tcga_wsRGES_cs_runs_used.tsv', sep='\t', index=False)
    out.to_csv(OUT_DIR / 'tcga_wsRGES_vs_IC50.tsv', sep='\t', index=False)

    # known answer: matched-line best / srges recomputed from the new runs vs the DEG rows of the registry
    k = res[(res.config == 'DEG') & res.collapse.isin(LINE_COLLAPSES)]
    k = k.merge(ref[ref.source == 'grid_DEG'], on=['cancer', 'collapse', 'ic50_cells', 'ic50_collapse'],
                suffixes=('', '_grid'))
    print('\nknown answer, DEG matched-line collapses (this run vs registry):')
    print(k[['cancer', 'collapse', 'ic50_cells', 'ic50_collapse', 'n_drugs', 'n_drugs_grid',
             'spearman_r', 'spearman_r_grid']].round(3).to_string(index=False))

    piv = out[out.config.isin(['DEG', 'Chen2017']) | out.source.str.startswith('grid')].copy()
    piv['method'] = np.where(piv.source == 'this_run', piv.collapse + '|' + piv.drug_set.astype(str),
                             piv.source + ':' + piv.config.astype(str) + ':' + piv.collapse.astype(str))
    piv = piv.pivot_table(index=['cancer', 'method'], columns=['ic50_cells', 'ic50_collapse'],
                          values='spearman_r', aggfunc='first')
    print('\nSpearman with IC50:')
    print(piv.round(2).to_string())
    boot = bootstrap_diffs(boot_in)
    boot.to_csv(OUT_DIR / 'tcga_wsRGES_bootstrap_spearman_diff.tsv', sep='\t', index=False)
    print('\npaired bootstrap of Spearman differences, matched line, median IC50 (DEG, line drugs):')
    print(boot[(boot.ic50_cells == 'matched') & (boot.ic50_collapse == 'median')]
          [['cancer', 'n_drugs', 'collapse_a', 'collapse_b', 'rho_a', 'rho_b', 'diff', 'ci_low', 'ci_high',
            'frac_boot_diff_le0']].round(3).to_string(index=False))
    figure(res, merged_fig, boot)
