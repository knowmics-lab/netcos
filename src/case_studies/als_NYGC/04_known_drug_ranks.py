#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Position of the known ALS drugs, and agreement between drug rankings, for the
drug-collapsed (TSR/LMM) connectivity-score runs and the per-instance sRGES
scores of 03_srges_per_drug.py, on the als_NYGC disease signature.

Design
------
- Rankings (one score per drug, more negative = stronger reversal):
  * '<config>_TSR': the most recent cs_batch run of the slug (= DISEASE id) on
    drug input 'all' (drugs collapsed before CS), per configuration in
    TSR_CONFIGS, resolved with loader.load_cs_run; score = pooled 6h_24h
    connectivity_score of the drug ('pooled' collapse) or the most negative of
    its 6h / 24h / 6h_24h scores ('best_timepoint' collapse).
  * '<config>_sRGES', '<config>_wsRGES', '<config>_best': the
    sRGES_unweighted, sRGES_weighted and best collapses of
    <slug>_srges_per_drug.tsv, for every configuration present there.
- Drug key: name in lower case, alphanumerics only (TSR: the drug part of
  pert_id; per-instance: pert_iname).
- Percentile of a drug = 100 x ascending rank / n within a ranking (low =
  reversal); 'own' over every drug of that ranking, 'shared' over the drugs
  scored by every DEG ranking (DEG_TSR, DEG_sRGES, DEG_wsRGES, DEG_best).
- Ranking agreement, for every pair of rankings (pooled collapse for TSR), on
  the drugs scored by both: Spearman correlation and overlap of the top-k
  (k in TOP_K) after re-ranking on the shared drugs.
- Known drugs: KNOWN_DRUGS (class labels only; no statistic uses them). For
  each, the number of L1000 instances, of instances on LINCS lines with a
  CCLE weight, and the lines without a weight that carry its instances.
- Lines without CCLE expression (no weight): their instances enter the
  unweighted sRGES and are left out of the weighted one; counted overall and
  per known drug, with the drugs scored by sRGES but not by wsRGES.
- The two routes use different LINCS instance pools: the '_TSR' runs come
  from GSE92742 level-3 profiles collapsed by a linear mixed model, the
  per-instance rankings from the level-5 gold set of the Chen 2017
  distribution; a comparison between them compares method and pool together.
- Instance-count strata: the per-instance scores of compounds with few
  instances are more extreme (and 'best' of compounds with many instances
  more negative), so for the per-instance rankings the known drugs also get a
  percentile within their stratum of instance count (STRATA, number of
  instances the ranking averages: 6h / 24h with dose > 0 for sRGES_unweighted
  and best, of those the ones on lines with a weight for sRGES_weighted),
  and the Spearman
  correlation of DEG_TSR with each DEG per-instance ranking is also given per
  stratum.
- Known-answer checks of the per-instance CS:
  (i) for N_CHECK random instances of CHECK_CELL_LINE (seed 0), the stored
  connectivity_score is recomputed with
  connectivity_score.bin_chen_connectivity on the same inputs;
  (ii) for every instance of every cell line, an independent implementation
  of the Chen 2017 RGES (cmap_score_new: KS statistics a, b of the up- and
  down-regulated disease genes on the positions of the drug ranking,
  ks = a if a > b else -b, score ks_up - ks_down, also when the two have the
  same sign) is computed with the 1-based positions of the published code
  and with the 0-based positions of connectivity_score.rank_genes (a + 1/r,
  b - 1/r), both on the drug ranking ordered as in rank_genes (descending
  pandas sort, so tied drug values keep the same order). The 0-based score is compared with the stored CS; the 1-based
  score is collapsed per compound with sRGES.compute_sRGES (unweighted and
  weighted) to measure the effect of the position convention on the
  per-compound ranking. The share of instances whose ks_up and ks_down have
  the same sign is reported (Lamb 2006 scores those instances 0; Chen 2017
  and this pipeline do not).

Inputs
------
- P.oout_dir/<slug>_srges_per_drug.tsv, <slug>_lincs_ccle_mapping.tsv
- conf.cs_log_filename and the CS files it points to
- conf.LINCS_METADATA_PATH
- case_conf.DISEASE_SIGNATURE and the per-instance DEG drug inputs
  (conf.CS_DIR/input/drug_signature_2025/<cell_id>_LM)
- P.oout_dir/<slug>_cell_line_weights.tsv

Outputs
-------
- P.oout_dir/<slug>_known_drug_ranks.tsv
    one row per (known drug, ranking, collapse): score, rank, n_drugs,
    percentile_own, percentile_shared, cs_run_id (TSR) and, for per-instance
    rankings, n_instances, n_cell_lines, instance_stratum,
    n_drugs_stratum, percentile_stratum; rankings '<...>_rank1based' are the
    sRGES / wsRGES recomputed from the 1-based per-instance scores of (ii).
- P.oout_dir/<slug>_known_drug_instances.tsv
    per known drug: L1000 instances (all, 6h / 24h with dose > 0, on lines
    with a weight), lines without a weight, in the TSR universe.
- P.oout_dir/<slug>_ranking_comparison.tsv
    one row per (pair of rankings, stratum): n_shared, spearman,
    top50_overlap, top100_overlap; stratum 'all' for every pair, instance-
    count strata for DEG_TSR vs the DEG per-instance rankings.
- P.oout_dir/<slug>_cs_known_answer_check.tsv
    per checked instance of (i): stored CS, CS recomputed with the pipeline
    function, CS of the independent implementation, ks_up, ks_down.
- P.oout_dir/<slug>_cs_position_convention.tsv
    per instance of (ii): LINCS_id, cell_id, cs_stored, cs_0based,
    cs_1based, ks_up_0based, ks_down_0based.
- P.img_dir/<slug>_tsr_vs_srges_panel.png
    A: percentile of the known drugs in every ranking; B: DEG_TSR vs
    DEG_sRGES percentile on the shared drugs, known drugs marked; C: Spearman
    correlation between rankings.

Usage
-----
    cd src && PYTHONPATH=. python case_studies/als_NYGC/04_known_drug_ranks.py
"""

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from case_conf import P, DISEASE_SIGNATURE
import conf
from loader import load_cs_run, load_single_signature_cs_input, standardize_disease_signature
from connectivity_score import get_common_genes, bin_chen_connectivity
from sRGES import compute_sRGES
from joblib import Parallel, delayed
import lfs_style

TSR_CONFIGS = {'DEG': (0, 0, 0), 'DEG_LM': (0, 1, 0), 'mith': (1, 0, 0), 'mith_LM': (1, 1, 0), 'mith_pw': (1, 0, 1)}
SRGES_NAMES = {'sRGES_unweighted': 'sRGES', 'sRGES_weighted': 'wsRGES', 'best': 'best'}
DEG_RANKINGS = ['DEG_TSR', 'DEG_sRGES', 'DEG_wsRGES', 'DEG_best']
KNOWN_DRUGS = {  # name: class
    'riluzole': 'approved, disease-modifying',
    'edaravone': 'approved, disease-modifying',
    'dextromethorphan': 'approved, symptomatic (Nuedexta component)',
    'quinidine': 'approved, symptomatic (Nuedexta component)',
    'phenylbutyrate': 'failed phase 3 (AMX0035 component)',
}
TIMEPOINTS = ('6h_24h', '24h', '6h')
TOP_K = (50, 100)
CHECK_CELL_LINE = 'NEU'
N_CHECK = 20
N_JOBS = 12
STRATA = ([0, 1, 2, 5, 9, 25, 50, np.inf], ['1', '2', '3-5', '6-9', '10-25', '26-50', '>50'])

lfs_style.apply(**{'savefig.dpi': 130, 'font.size': 10, 'legend.shadow': False})


def key(s):
    return pd.Series(s, dtype=str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True).values


def load_tsr():
    """{'<config>_TSR': DataFrame(drug_key, pooled, best_timepoint)}, {label: cs_run_id}."""
    out, ids = {}, {}
    for c, (m, lm, pw) in TSR_CONFIGS.items():
        df, path, run_id = load_cs_run(disease=P.slug, landmark_disease=False, drug_run_id='all', mith=m,
                                       cs_on_LM=lm, CS_ON_PATHWAYS=pw, CS_METHOD='bin_chen', verbose=False)
        tp = df['pert_id'].str.extract(rf'^(.*?)_({"|".join(TIMEPOINTS)})$')
        df['drug_key'] = key(tp[0])
        df['tp'] = tp[1]
        pooled = df[df.tp == '6h_24h'].groupby('drug_key')['connectivity_score'].mean()
        best = df.groupby('drug_key')['connectivity_score'].min()
        out[f'{c}_TSR'] = pd.DataFrame({'pooled': pooled, 'best_timepoint': best})
        ids[f'{c}_TSR'] = run_id
        print(f'{c}_TSR: {run_id}; {df.drug_key.nunique()} drugs, {len(df)} rows')
    return out, ids


def load_per_instance():
    t = pd.read_csv(P.oout_dir / f'{P.slug}_srges_per_drug.tsv', sep='\t')
    t['drug_key'] = key(t['pert_iname'])
    out = {}
    for (c, col), g in t.groupby(['config', 'collapse']):
        g = g.groupby('drug_key').agg(score=('score', 'mean'), n_instances=('n_instances', 'sum'),
                                       n_cell_lines=('n_cell_lines', 'max'))
        out[f'{c}_{SRGES_NAMES[col]}'] = g
    return out


def drug_dir(cell):
    return conf.CS_DIR / 'input' / 'drug_signature_2025' / (cell + '_LM')


def disease_common():
    """DEG disease signature restricted to the genes of the drug signatures."""
    dis = pd.read_csv(DISEASE_SIGNATURE, sep=';', decimal=',', dtype={'gene_id': str})
    dis = standardize_disease_signature(dis, source='tsr')
    first = next(drug_dir(CHECK_CELL_LINE).glob('*_signature_gene_id.pkl'))
    drug = pd.read_pickle(first)
    i_dis, _ = get_common_genes(dis[['gene_id', 'DE_log2_FC', 'adj.p.value']], drug[['gene_id', 'DE_log2_FC', 'adj.p.value']])
    return dis.iloc[i_dis].reset_index(drop=True)


def _ab(v, r):
    s = len(v)
    j = np.arange(1, s + 1)
    return np.max(j / s - v / r), np.max(v / r - (j - 1) / s)


def _ks(a, b):
    return a if a > b else -b


def ks_one_line(cell, ids_cs, up, down):
    rows = []
    for i, cs in ids_cs:
        d = load_single_signature_cs_input(i, drug_dir(cell), mith=False)
        r = len(d)
        # drug genes in the order of connectivity_score.rank_genes (descending pandas sort, same tie order)
        order = d.sort_values(by='DE_log2_FC', ascending=False)['gene_id'].values
        pos = pd.Series(np.arange(1, r + 1), index=order)
        (au, bu), (ad, bd) = (_ab(np.sort(pos.loc[g].values), r) for g in (up, down))
        k0u, k0d = _ks(au + 1 / r, bu - 1 / r), _ks(ad + 1 / r, bd - 1 / r)
        rows.append({'LINCS_id': i, 'cell_id': cell, 'cs_stored': cs, 'cs_0based': k0u - k0d,
                     'cs_1based': _ks(au, bu) - _ks(ad, bd), 'ks_up_0based': k0u, 'ks_down_0based': k0d})
    return rows


def position_convention_check(md, weights):
    dis = disease_common()
    up = dis.loc[dis.DE_log2_FC > 0, 'gene_id'].values
    down = dis.loc[dis.DE_log2_FC < 0, 'gene_id'].values
    jobs = []
    for c in sorted(md.cell_id.unique()):
        df, _, _ = load_cs_run(disease=P.slug, landmark_disease=False, drug_run_id=c + '_LM', mith=0, cs_on_LM=0,
                               CS_ON_PATHWAYS=0, CS_METHOD='bin_chen', verbose=False)
        jobs.append((c, list(zip(df['LINCS_id'].astype(str), df['connectivity_score']))))
    res = Parallel(n_jobs=N_JOBS)(delayed(ks_one_line)(c, ids, up, down) for c, ids in jobs)
    ch = pd.DataFrame([r for rr in res for r in rr])
    d0 = (ch.cs_stored - ch.cs_0based).abs()
    d1 = (ch.cs_1based - ch.cs_0based).abs()
    print(f'\nposition convention, all {len(ch)} instances ({len(up)} up, {len(down)} down disease genes):')
    print(f'  max |stored - independent 0-based| = {d0.max():.2e}')
    print(f'  1-based vs 0-based: differ in {(d1 > 1e-12).mean():.1%} of instances, by > 0.01 in '
          f'{(d1 > 0.01).mean():.1%}; max {d1.max():.4f}')
    print(f'  instances with ks_up and ks_down of the same sign: {(np.sign(ch.ks_up_0based) == np.sign(ch.ks_down_0based)).mean():.1%}')
    x = ch.drop(columns='cell_id').merge(md, left_on='LINCS_id', right_on='id')
    out = {}
    for name, w in (('sRGES', None), ('wsRGES', weights)):
        a = compute_sRGES(x, score_col='cs_stored', drug_col='pert_iname', cell_col='cell_id',
                          cell_line_weights=w, srges_diff_mode='binchen_r').set_index('pert_iname')['sRGES']
        b = compute_sRGES(x, score_col='cs_1based', drug_col='pert_iname', cell_col='cell_id',
                          cell_line_weights=w, srges_diff_mode='binchen_r').set_index('pert_iname')['sRGES']
        dd = pd.concat([a, b], axis=1, keys=['stored', 'r1']).dropna()
        print(f'  DEG_{name}, stored vs 1-based per-instance CS: Spearman {spearmanr(dd.stored, dd.r1).statistic:.4f}, '
              f'max |diff| {(dd.stored - dd.r1).abs().max():.4f}')
        b.index = key(b.index)
        out[f'DEG_{name}_rank1based'] = b.groupby(level=0).mean()
    return ch, out


def percentile(s):
    return s.rank(method='average') / s.notna().sum() * 100


def known_answer_check(md):
    m = md[md.cell_id == CHECK_CELL_LINE]
    df, path, run_id = load_cs_run(disease=P.slug, landmark_disease=False, drug_run_id=CHECK_CELL_LINE + '_LM',
                                   mith=0, cs_on_LM=0, CS_ON_PATHWAYS=0, CS_METHOD='bin_chen', verbose=False)
    df['LINCS_id'] = df['LINCS_id'].astype(str)
    ids = df['LINCS_id'].sample(n=min(N_CHECK, len(df)), random_state=0).tolist()
    dis = pd.read_csv(DISEASE_SIGNATURE, sep=';', decimal=',', dtype={'gene_id': str})
    dis = standardize_disease_signature(dis, source='tsr')
    ddir = drug_dir(CHECK_CELL_LINE)
    cols = ['gene_id', 'DE_log2_FC', 'adj.p.value']
    first = load_single_signature_cs_input(ids[0], ddir, mith=False)
    i_dis, _ = get_common_genes(dis[cols], first[cols])
    dis_c = dis[cols].iloc[i_dis].reset_index(drop=True)
    rows = []
    for i in ids:
        drug = load_single_signature_cs_input(i, ddir, mith=False)[cols]
        cs_pipe, _ = bin_chen_connectivity(dis_c, drug)
        # independent RGES: 1-based positions of the disease genes in the drug ranking (most up = 1)
        pos = pd.Series(drug['DE_log2_FC'].rank(ascending=False, method='first').values, index=drug['gene_id'].values)
        r = len(drug)
        ks = {}
        for side, sel in (('up', dis_c['DE_log2_FC'] > 0), ('down', dis_c['DE_log2_FC'] < 0)):
            v = np.sort(pos.loc[dis_c.loc[sel, 'gene_id']].values)
            s = len(v)
            a = np.max(np.arange(1, s + 1) / s - v / r)
            b = np.max(v / r - np.arange(s) / s)
            ks[side] = a if a > b else -b
        rows.append({'cs_run_id': run_id, 'LINCS_id': i,
                     'pert_iname': m.set_index('id').loc[i, 'pert_iname'] if i in set(m['id']) else None,
                     'cs_stored': df.set_index('LINCS_id').loc[i, 'connectivity_score'],
                     'cs_pipeline_recomputed': cs_pipe, 'cs_independent': ks['up'] - ks['down'],
                     'ks_up': ks['up'], 'ks_down': ks['down']})
    ch = pd.DataFrame(rows)
    d1 = (ch.cs_stored - ch.cs_pipeline_recomputed).abs().max()
    d2 = (ch.cs_stored - ch.cs_independent).abs().max()
    same = (np.sign(ch.ks_up) == np.sign(ch.ks_down)).mean()
    print(f'\nknown-answer check, {CHECK_CELL_LINE}_LM ({run_id}), {len(ch)} instances, '
          f'{len(dis_c)} disease genes on the drug genes ({(dis_c.DE_log2_FC > 0).sum()} up, '
          f'{(dis_c.DE_log2_FC < 0).sum()} down), r = {r}:')
    print(f'  max |stored - pipeline recomputed| = {d1:.2e}')
    print(f'  max |stored - independent RGES|    = {d2:.2e}  (1/r = {1 / r:.2e}); '
          f'Spearman {spearmanr(ch.cs_stored, ch.cs_independent).statistic:.4f}')
    print(f'  instances with ks_up and ks_down of the same sign: {same:.0%}')
    return ch


if __name__ == '__main__':
    md = pd.read_csv(conf.LINCS_METADATA_PATH, dtype={'id': str})
    mapping = pd.read_csv(P.oout_dir / f'{P.slug}_lincs_ccle_mapping.tsv', sep='\t')
    no_w = set(mapping.loc[~mapping.has_expression, 'lincs_cell_id'])

    w = pd.read_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t')
    weights = dict(zip(w['lincs_cell_id'], w['cor']))
    tsr, tsr_ids = load_tsr()
    inst = load_per_instance()
    # instance-count stratum of every compound in each per-instance ranking: the instances that ranking
    # averages (sRGES_unweighted and best: all 6h / 24h instances with dose > 0; sRGES_weighted: those on
    # lines with a weight)
    strata = {k: pd.cut(v['n_instances'], STRATA[0], labels=STRATA[1]).astype(str) for k, v in inst.items()}
    print('per-instance rankings:', ', '.join(f'{k} ({len(v)})' for k, v in inst.items()))

    # name matching between the two routes
    md['drug_key'] = key(md['pert_iname'])
    un = set(tsr['DEG_TSR'].index) - set(md['drug_key'])
    print(f'DEG_TSR drugs whose name is not an L1000 pert_iname of the metadata: {len(un)}'
          + (f' (e.g. {", ".join(sorted(un)[:8])})' if un else ''))

    # rankings as score series (pooled collapse for TSR)
    R = {k: v['pooled'].dropna() for k, v in tsr.items()}
    R.update({k: v['score'].dropna() for k, v in inst.items()})
    order = [k for k in [f'{c}_{s}' for c in TSR_CONFIGS for s in ('TSR', 'sRGES', 'wsRGES', 'best')] if k in R]
    shared = set.intersection(*[set(R[k].index) for k in DEG_RANKINGS if k in R])
    print('\ndrug universes: ' + ', '.join(f'{k} {len(R[k])}' for k in order) +
          f'; shared by the DEG rankings {len(shared)}')

    # instances on lines without a weight
    kept = md[(pd.to_numeric(md.pert_dose, errors='coerce') > 0) & md.pert_time.isin([6, 24])]
    n_nw = kept.cell_id.isin(no_w).sum()
    print(f'\n6h/24h instances with dose > 0: {len(kept)}; on the {len(no_w)} lines without a weight: {n_nw} '
          f'({n_nw / len(kept):.1%})')
    if 'DEG_sRGES' in R and 'DEG_wsRGES' in R:
        lost = set(R['DEG_sRGES'].index) - set(R['DEG_wsRGES'].index)
        print(f'drugs scored by DEG_sRGES but not by DEG_wsRGES (every instance on a line without weight): {len(lost)}')

    # known drugs
    krows, irows = [], []
    for d, cls_ in KNOWN_DRUGS.items():
        k = key([d])[0]
        mi = md[md.drug_key == k]
        mk = kept[kept.drug_key == k]
        irows.append({'drug': d, 'class': cls_, 'n_instances_all': len(mi), 'n_instances_6h24h': len(mk),
                      'n_instances_weighted_lines': int((~mk.cell_id.isin(no_w)).sum()),
                      'times_h': ','.join(map(str, sorted(mk.pert_time.unique()))),
                      'lines_without_weight': ','.join(sorted(set(mk.cell_id) & no_w)),
                      'lines_with_weight': ','.join(sorted(set(mk.cell_id) - no_w)),
                      'in_TSR_universe': k in tsr['DEG_TSR'].index})
        for name in order:
            colls = ['pooled', 'best_timepoint'] if name.endswith('_TSR') else [None]
            for coll in colls:
                s = tsr[name][coll].dropna() if coll else R[name]
                rec = {'drug': d, 'class': cls_, 'ranking': name, 'collapse': coll or name.split('_')[-1],
                       'n_drugs': len(s), 'cs_run_id': tsr_ids.get(name)}
                if k in s.index:
                    rec.update(score=s[k], rank=int(s.rank(method='min')[k]), percentile_own=percentile(s)[k],
                               percentile_shared=percentile(s[s.index.isin(shared)])[k] if k in shared else np.nan)
                    if name in inst:
                        rec.update(n_instances=inst[name].loc[k, 'n_instances'], n_cell_lines=inst[name].loc[k, 'n_cell_lines'])
                        stratum = strata[name]
                        st = s[stratum.reindex(s.index) == stratum[k]]
                        rec.update(instance_stratum=stratum[k], n_drugs_stratum=len(st), percentile_stratum=percentile(st)[k])
                krows.append(rec)
    kn = pd.DataFrame(krows)
    ki = pd.DataFrame(irows)
    print('\nknown drugs, L1000 instances:')
    print(ki.to_string(index=False))
    pv = kn[kn.collapse.isin(['TSR', 'sRGES', 'wsRGES', 'best']) | (kn.collapse == 'pooled')]
    print('\nknown drugs, percentile_own (low = reversal; TSR = pooled 6h_24h):')
    print(pv.pivot_table(index='drug', columns='ranking', values='percentile_own', dropna=False).reindex(columns=order).round(1).to_string())
    print('\nknown drugs, percentile on the drugs shared by the DEG rankings:')
    print(pv.pivot_table(index='drug', columns='ranking', values='percentile_shared', dropna=False)
          .reindex(columns=[c for c in DEG_RANKINGS if c in R]).round(1).to_string())
    print('\nknown drugs, TSR best-timepoint percentile_own:')
    print(kn[kn.collapse == 'best_timepoint'].pivot_table(index='drug', columns='ranking', values='percentile_own')
          .round(1).to_string())

    # ranking agreement
    crow = []
    for a, b in combinations(order, 2):
        s = pd.concat([R[a], R[b]], axis=1, join='inner')
        rec = {'ranking_a': a, 'ranking_b': b, 'n_shared': len(s),
               'spearman': spearmanr(s.iloc[:, 0], s.iloc[:, 1]).statistic}
        for kk in TOP_K:
            rec[f'top{kk}_overlap'] = len(set(s.iloc[:, 0].nsmallest(kk).index) & set(s.iloc[:, 1].nsmallest(kk).index))
        crow.append(rec)
    for b in [c for c in DEG_RANKINGS[1:] if c in R]:
        s = pd.concat([R['DEG_TSR'], R[b]], axis=1, join='inner')
        for st in STRATA[1]:
            ss = s[strata[b].reindex(s.index) == st]
            if len(ss) < 3:
                continue
            crow.append({'ranking_a': 'DEG_TSR', 'ranking_b': b, 'stratum': st, 'n_shared': len(ss),
                         'spearman': spearmanr(ss.iloc[:, 0], ss.iloc[:, 1]).statistic})
    comp = pd.DataFrame(crow)
    comp['stratum'] = comp['stratum'].fillna('all')
    print('\nranking agreement (shared drugs):')
    print(comp[comp.stratum == 'all'].round(3).to_string(index=False))
    print('\nDEG_TSR vs DEG per-instance rankings, by instance-count stratum (instances averaged by the per-instance ranking):')
    print(comp[comp.stratum != 'all'].pivot_table(index='stratum', columns='ranking_b', values='spearman')
          .reindex(STRATA[1]).round(3).to_string())
    print(comp[comp.stratum != 'all'].pivot_table(index='stratum', columns='ranking_b', values='n_shared')
          .reindex(STRATA[1]).to_string())
    print('\nknown drugs, percentile within the instance-count stratum (per-instance rankings):')
    print(kn[kn.percentile_stratum.notna()].pivot_table(index=['drug', 'instance_stratum', 'n_drugs_stratum'],
          columns='ranking', values='percentile_stratum').round(1).to_string())

    ch = known_answer_check(md)
    conv, r1 = position_convention_check(md, weights)
    extra = []
    for name, s in r1.items():
        for d, cls_ in KNOWN_DRUGS.items():
            k = key([d])[0]
            if k in s.index:
                extra.append({'drug': d, 'class': cls_, 'ranking': name, 'collapse': name.split('_')[1],
                              'n_drugs': len(s), 'score': s[k], 'rank': int(s.rank(method='min')[k]),
                              'percentile_own': percentile(s)[k]})
    kn = pd.concat([kn, pd.DataFrame(extra)], ignore_index=True)
    print('\nknown drugs, percentile_own with the 1-based per-instance CS:')
    print(pd.DataFrame(extra).pivot_table(index='drug', columns='ranking', values='percentile_own').round(1).to_string())

    P.oout_dir.mkdir(parents=True, exist_ok=True)
    kn.to_csv(P.oout_dir / f'{P.slug}_known_drug_ranks.tsv', sep='\t', index=False)
    ki.to_csv(P.oout_dir / f'{P.slug}_known_drug_instances.tsv', sep='\t', index=False)
    comp.to_csv(P.oout_dir / f'{P.slug}_ranking_comparison.tsv', sep='\t', index=False)
    ch.to_csv(P.oout_dir / f'{P.slug}_cs_known_answer_check.tsv', sep='\t', index=False)
    conv.to_csv(P.oout_dir / f'{P.slug}_cs_position_convention.tsv', sep='\t', index=False)

    # panel
    fig = plt.figure(figsize=(16, 11))
    gsp = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], wspace=0.35)
    a = fig.add_subplot(gsp[0, :])
    cols = plt.cm.tab10.colors
    for j, name in enumerate(order):
        sub = pv[pv.ranking == name].set_index('drug')
        y = np.arange(len(KNOWN_DRUGS)) + (j - len(order) / 2) * 0.045
        x = [sub.loc[d, 'percentile_own'] if d in sub.index else np.nan for d in KNOWN_DRUGS]
        mk = 's' if name.endswith('_TSR') else ('o' if name.endswith('sRGES') else '^')
        a.scatter(x, y, marker=mk, s=30, color=cols[j % 10], label=name,
                  edgecolor='k' if name.startswith('DEG') else 'none', linewidth=.5)
    a.set_yticks(np.arange(len(KNOWN_DRUGS)))
    a.set_yticklabels([f'{d}\n({c.split(",")[0]})' for d, c in KNOWN_DRUGS.items()], fontsize=9)
    a.set_xlim(0, 100)
    a.invert_yaxis()
    a.set_xlabel('percentile within the ranking (low = stronger reversal)')
    a.set_title('A. Known ALS drugs: percentile in each ranking (squares _TSR, circles sRGES/wsRGES, '
                'triangles best; missing = not scored)')
    a.legend(fontsize=8, ncol=2, loc='center left', bbox_to_anchor=(1.0, .5))

    b = fig.add_subplot(gsp[1, 0])
    if 'DEG_sRGES' in R:
        s = pd.concat([R['DEG_TSR'], R['DEG_sRGES']], axis=1, join='inner', keys=['t', 's'])
        pt, ps = percentile(s.t), percentile(s.s)
        b.scatter(pt, ps, s=4, color='0.6', alpha=.5)
        for d in KNOWN_DRUGS:
            k = key([d])[0]
            if k in s.index:
                b.scatter(pt[k], ps[k], s=40, color=lfs_style.AFTER)
                b.annotate(d, (pt[k], ps[k]), xytext=(4, 4), textcoords='offset points', fontsize=9)
        rho = spearmanr(s.t, s.s).statistic
        b.set_xlabel('DEG_TSR percentile (pooled 6h_24h CS)')
        b.set_ylabel('DEG_sRGES percentile (unweighted sRGES)')
        b.set_title(f'B. DEG_TSR vs DEG_sRGES on {len(s)} shared drugs (Spearman {rho:.2f})')

    c = fig.add_subplot(gsp[1, 1])
    M = pd.DataFrame(np.eye(len(order)), index=order, columns=order)
    for _, r in comp.iterrows():
        M.loc[r.ranking_a, r.ranking_b] = M.loc[r.ranking_b, r.ranking_a] = r.spearman
    im = c.imshow(M.values, vmin=-1, vmax=1, cmap='RdBu_r')
    c.set_xticks(range(len(order)))
    c.set_xticklabels(order, rotation=60, ha='right', fontsize=8)
    c.set_yticks(range(len(order)))
    c.set_yticklabels(order, fontsize=8)
    for i in range(len(order)):
        for j in range(len(order)):
            c.text(j, i, f'{M.values[i, j]:.2f}', ha='center', va='center', fontsize=7)
    c.grid(False)
    c.set_title('C. Spearman correlation between rankings (shared drugs)')
    fig.colorbar(im, ax=c, fraction=.046)
    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_tsr_vs_srges_panel.png')
