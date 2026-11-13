#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Drug rankings of ipf from the drug-collapsed (TSR/LMM) connectivity-score runs
and from the per-instance sRGES scores of 03_srges_per_drug.py: position of
the known IPF drugs, agreement between rankings, size of the drug universes.

Design
------
- Rankings (one score per drug, more negative = stronger reversal):
  * '<config>_TSR': the pooled 6h_24h connectivity_score of the most recent
    cs_batch run of conf.DISEASE on drug input 'all' (drugs collapsed by the
    TSR/LMM model before CS), per configuration in TSR_CONFIGS, resolved with
    loader.load_cs_run.
  * '<config>_sRGES', '<config>_wsRGES', '<config>_best': the
    'sRGES_unweighted', 'sRGES_weighted' and 'best' per-drug scores of
    03_srges_per_drug.py, for every configuration present in its output.
- Drug key: the L1000 compound name in lower case, alphanumerics only
  (pert_iname in the per-instance route, the drug part of pert_id in the TSR
  route); names that collapse to the same key take the mean score.
- Percentile of a drug in a ranking: rank(pct=True) x 100 of its score in
  ascending order (low percentile = strong reverser), over the full universe
  of that ranking and over the drugs scored by both DEG_TSR and DEG_sRGES
  ('shared_DEG'; for a ranking that scores only part of them, e.g.
  DEG_wsRGES, over that part).
- Instance-count strata (per-instance rankings): the drug's number of
  instances binned by STRATA; 'percentile_stratum' is the percentile among
  the drugs of the same ranking and stratum (scores of drugs with few
  instances are more dispersed, and 'best' falls with the number of
  instances).
- Instance pools: the TSR route collapses the LINCS GSE92742 level-3
  instances of a drug (per-drug LMM of tsr/); the per-instance route uses the
  level-5 gold set of conf.LINCS_METADATA_PATH. A TSR vs per-instance
  comparison therefore differs in both the collapse and the instances.
- Known IPF drugs: KNOWN_DRUGS (approved, disease-modifying). For the
  per-instance route the instance-level CS of these drugs (cell line, dose,
  time) is read from the runs listed in <slug>_srges_cs_runs_used.tsv.
- Cross-check of the TSR route: the same percentiles with the drug score
  taken as the mean of the 6h, 24h and 6h_24h rows ('TSR_mean3tp').
- Agreement of two rankings on their shared drugs: Spearman correlation and
  overlap of the top-k (k in TOP_K), the top-k taken after restricting both
  rankings to the shared drugs; for every TSR ranking vs DEG_sRGES the same
  within each instance-count stratum of DEG_sRGES (column 'stratum';
  'all' = no stratification).
- Universe sizes: drugs per ranking; per-instance drugs absent from the TSR
  universe, and how many of those have only 6 h instances (dose > 0, 6 or
  24 h, the instances sRGES keeps).

conf.DISEASE must equal the slug (the script stops otherwise).

Inputs
------
- P.oout_dir/<slug>_srges_per_drug.tsv, P.oout_dir/<slug>_srges_cs_runs_used.tsv
- conf.cs_log_filename and the CS files it points to
- conf.LINCS_METADATA_PATH

Outputs
-------
- P.oout_dir/<slug>_drug_rankings.tsv
    one row per (ranking, drug_key): drug, score, rank, percentile,
    percentile_shared_DEG, n_instances, n_cell_lines, stratum,
    percentile_stratum (per-instance only), cs_run_id (TSR only).
- P.oout_dir/<slug>_known_drug_ranks.tsv
    one row per (known drug, ranking): score, rank, n_drugs, percentile,
    percentile_shared_DEG, n_instances, n_cell_lines, stratum,
    n_drugs_stratum, percentile_stratum; rows 'TSR_mean3tp' for the
    cross-check.
- P.oout_dir/<slug>_ranking_agreement.tsv
    one row per pair of rankings (and stratum): stratum, n_a, n_b,
    n_shared, spearman, top<k>_overlap for every k in TOP_K.
- P.oout_dir/<slug>_known_drug_instances.tsv
    the per-instance CS of the known drugs (config, cell_id, pert_id,
    pert_dose, pert_time, connectivity_score, cell-line weight).
- P.img_dir/<slug>_srges_vs_tsr_panel.png
    A: instances and cell lines per compound (per-instance route, DEG);
    B: score distributions of the DEG rankings; C: DEG_TSR vs DEG_sRGES on
    the shared drugs; D: DEG_wsRGES vs DEG_sRGES; E: percentile of the known
    drugs in every ranking (own universe, shared drugs, instance-count
    stratum); F: Spearman correlation between rankings on the
    shared drugs.

Usage
-----
    cd src && PYTHONPATH=. python case_studies/ipf/04_known_drug_ranks.py
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

from case_conf import P
import conf
from loader import load_cs_run
import lfs_style

KNOWN_DRUGS = ['pirfenidone', 'nintedanib']
TSR_CONFIGS = {'DEG': (0, 0, 0), 'DEG_LM': (0, 1, 0), 'mith': (1, 0, 0), 'mith_LM': (1, 1, 0), 'mith_pw': (1, 0, 1)}
PER_INSTANCE = {'sRGES_unweighted': 'sRGES', 'sRGES_weighted': 'wsRGES', 'best': 'best'}
TOP_K = [50, 100]
STRATA = [0, 1, 2, 5, 9, 25, 50, np.inf]
STRATA_LABELS = ['1', '2', '3-5', '6-9', '10-25', '26-50', '>50']
TIMEPOINTS = ('6h_24h', '24h', '6h')
BLUE, RED = lfs_style.BEFORE, lfs_style.AFTER


def key(s):
    return pd.Series(s).astype(str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True).values


def split_pert_id(pid):
    for tp in TIMEPOINTS:
        if pid.endswith('_' + tp):
            return pid[:-(len(tp) + 1)], tp
    return pid, None


def tsr_rankings():
    out, mean3 = {}, {}
    for label, (m, lm, pw) in TSR_CONFIGS.items():
        try:
            df, path, run_id = load_cs_run(disease=conf.DISEASE, landmark_disease=False, drug_run_id='all', mith=m,
                                           cs_on_LM=lm, CS_ON_PATHWAYS=pw, CS_METHOD=conf.CS_METHOD, verbose=False)
        except FileNotFoundError:
            print(f'{label}_TSR: no run')
            continue
        dt = df['pert_id'].astype(str).map(split_pert_id)
        df = df.assign(drug=dt.str[0], tp=dt.str[1])
        df['drug_key'] = key(df['drug'])
        pooled = df[df['tp'] == '6h_24h']
        r = pooled.groupby('drug_key').agg(score=('connectivity_score', 'mean'), drug=('drug', 'first')).reset_index()
        r['cs_run_id'] = run_id
        out[f'{label}_TSR'] = r
        g = df.groupby('drug_key').agg(score=('connectivity_score', 'mean'), drug=('drug', 'first')).reset_index()
        g['cs_run_id'] = run_id
        mean3[f'{label}_TSR_mean3tp'] = g
        print(f'{label}_TSR: {run_id}; {len(r)} drugs (6h_24h), {len(g)} drugs (any timepoint)')
    return out, mean3


def per_instance_rankings():
    s = pd.read_csv(P.oout_dir / f'{P.slug}_srges_per_drug.tsv', sep='\t')
    s['drug_key'] = key(s['pert_iname'])
    out = {}
    for (config, coll), d in s.groupby(['config', 'collapse']):
        r = d.groupby('drug_key').agg(score=('score', 'mean'), drug=('pert_iname', 'first'),
                                      n_instances=('n_instances', 'sum'),
                                      n_cell_lines=('n_cell_lines', 'max')).reset_index()
        out[f'{config}_{PER_INSTANCE[coll]}'] = r
    return out


def add_ranks(r):
    r = r.copy()
    r['rank'] = r['score'].rank(method='min')
    r['percentile'] = r['score'].rank(pct=True) * 100
    return r


def agreement(rk, names, stratum='all', keep=None):
    rows = []
    for a, b in combinations(names, 2):
        x = rk[a].set_index('drug_key')['score']
        y = rk[b].set_index('drug_key')['score']
        sh = x.index.intersection(y.index)
        if keep is not None:
            sh = sh.intersection(keep)
        row = {'ranking_a': a, 'ranking_b': b, 'stratum': stratum, 'n_a': len(x), 'n_b': len(y), 'n_shared': len(sh),
               'spearman': spearmanr(x[sh], y[sh]).statistic if len(sh) > 2 else np.nan}
        for k in TOP_K:
            row[f'top{k}_overlap'] = len(set(x[sh].nsmallest(k).index) & set(y[sh].nsmallest(k).index))
        rows.append(row)
    return pd.DataFrame(rows)


def known_instances(md, weights):
    used = pd.read_csv(P.oout_dir / f'{P.slug}_srges_cs_runs_used.tsv', sep='\t').dropna(subset=['cs_run_id'])
    ids = md[pd.Series(key(md['pert_iname']), index=md.index).isin(KNOWN_DRUGS)]
    rows = []
    for (config, cell), run_id in used.set_index(['config', 'cell_line'])['cs_run_id'].items():
        sub = ids[ids['cell_id'] == cell]
        if sub.empty:
            continue
        df, _, _ = load_cs_run(selected_cs_run_id=run_id, disease=conf.DISEASE, landmark_disease=False, verbose=False)
        df = df.assign(LINCS_id=df['LINCS_id'].astype(str))
        x = sub.merge(df[['LINCS_id', 'connectivity_score']], left_on='id', right_on='LINCS_id')
        x.insert(0, 'config', config)
        rows.append(x)
    x = pd.concat(rows, ignore_index=True)
    x['weight'] = x['cell_id'].map(weights)
    return x[['config', 'pert_iname', 'cell_id', 'pert_id', 'id', 'pert_dose', 'pert_time',
              'connectivity_score', 'weight']].sort_values(['config', 'pert_iname', 'connectivity_score'])


def panel(rk, ag, kn, deg_names):
    fig, ax = plt.subplots(2, 3, figsize=(19, 11.5))
    s = rk['DEG_sRGES']

    a = ax[0, 0]
    bins = np.logspace(0, np.log10(max(s['n_instances'].max(), 2)) + .05, 30)
    a.hist(s['n_instances'], bins=bins, histtype='step', lw=1.6, color=BLUE, label='instances per compound')
    a.hist(s['n_cell_lines'], bins=bins, histtype='step', lw=1.6, color=RED, label='cell lines per compound')
    a.set_xscale('log'); a.set_yscale('log')
    a.set_xlabel('count'); a.set_ylabel('compounds')
    a.set_title(f'A. Per-instance route (DEG): {len(s):,} compounds')
    a.legend()

    a = ax[0, 1]
    cols = {'DEG_sRGES': BLUE, 'DEG_wsRGES': 'k', 'DEG_best': '0.5', 'DEG_TSR': RED}
    lo = min(rk[n]['score'].quantile(.001) for n in cols)
    hi = max(rk[n]['score'].quantile(.999) for n in cols)
    bins = np.linspace(lo, hi, 120)
    for n, c in cols.items():
        v = rk[n]['score']
        a.hist(v, bins=bins, histtype='step', lw=1.4, color=c, density=True,
               label=f'{n} (n={len(v):,}, median {v.median():+.3f}, {100 * (v < 0).mean():.0f}% < 0)')
    a.axvline(0, color='k', lw=.7)
    a.set_xlabel('per-drug score'); a.set_ylabel('density')
    a.set_title('B. Score distributions (DEG)')
    a.legend(fontsize=8)

    def scatter(a, nx, ny, title):
        x = rk[nx].set_index('drug_key')['score']
        y = rk[ny].set_index('drug_key')['score']
        sh = x.index.intersection(y.index)
        rho = spearmanr(x[sh], y[sh]).statistic
        a.scatter(x[sh], y[sh], s=3, color='0.55', alpha=.4, rasterized=True)
        for d, c in zip(KNOWN_DRUGS, [BLUE, RED]):
            if d in sh:
                a.scatter(x[d], y[d], s=60, color=c, edgecolor='k', zorder=3, label=d)
        a.axhline(0, color='k', lw=.5); a.axvline(0, color='k', lw=.5)
        a.set_xlabel(nx); a.set_ylabel(ny)
        a.set_title(f'{title} (n={len(sh):,}, Spearman {rho:.2f})')
        a.legend(fontsize=8, loc='upper left')

    scatter(ax[0, 2], 'DEG_TSR', 'DEG_sRGES', 'C. DEG_TSR vs DEG_sRGES')
    scatter(ax[1, 0], 'DEG_sRGES', 'DEG_wsRGES', 'D. DEG_sRGES vs DEG_wsRGES')

    a = ax[1, 1]
    order = [n for n in rk if n.endswith('_TSR')] + [n for n in rk if not n.endswith('_TSR')]
    xpos = np.arange(len(order))
    for d, c, dx in zip(KNOWN_DRUGS, [BLUE, RED], [-.12, .12]):
        k = kn[kn['drug_key'] == d].set_index('ranking')
        a.scatter(xpos + dx, k['percentile'].reindex(order), color=c, s=45, zorder=3, label=f'{d}, own universe')
        a.scatter(xpos + dx, k['percentile_shared_DEG'].reindex(order), facecolor='none', edgecolor=c, s=45,
                  zorder=3, label=f'{d}, DEG_TSR & DEG_sRGES drugs')
        a.scatter(xpos + dx, k['percentile_stratum'].reindex(order), color=c, marker='^', s=35, alpha=.6,
                  zorder=3, label=f'{d}, same instance-count stratum')
    a.axhline(50, color='k', lw=.6)
    a.set_xticks(xpos); a.set_xticklabels(order, rotation=45, ha='right')
    a.set_ylim(0, 100); a.set_ylabel('percentile (low = stronger reversal)')
    a.set_title('E. Known IPF drugs')
    a.legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=2)

    a = ax[1, 2]
    names = deg_names + [n for n in rk if n.endswith('_TSR') and n not in deg_names]
    M = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for _, r in ag.iterrows():
        if r.stratum == 'all' and r.ranking_a in names and r.ranking_b in names:
            M.loc[r.ranking_a, r.ranking_b] = M.loc[r.ranking_b, r.ranking_a] = r.spearman
    im = a.imshow(M.values, vmin=-1, vmax=1, cmap='RdBu_r')
    a.set_xticks(range(len(names))); a.set_xticklabels(names, rotation=45, ha='right')
    a.set_yticks(range(len(names))); a.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            a.text(j, i, f'{M.values[i, j]:.2f}', ha='center', va='center', fontsize=8)
    a.grid(False)
    a.set_title('F. Spearman between rankings (shared drugs)')
    fig.colorbar(im, ax=a, fraction=.046)

    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_srges_vs_tsr_panel.png', dpi=130)


if __name__ == '__main__':
    assert conf.DISEASE == P.slug, f'conf.DISEASE is {conf.DISEASE!r}; this script needs the ipf conf'
    lfs_style.apply()
    md = pd.read_csv(conf.LINCS_METADATA_PATH, dtype={'id': str})
    md = md[['id', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose', 'pert_time']]
    w = pd.read_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t')
    weights = dict(zip(w['lincs_cell_id'], w['cor']))

    tsr, mean3 = tsr_rankings()
    pi = per_instance_rankings()
    rk = {n: add_ranks(r) for n, r in {**tsr, **pi}.items()}
    mean3 = {n: add_ranks(r) for n, r in mean3.items()}
    deg_names = [n for n in rk if n.startswith('DEG_') and not n.startswith('DEG_LM')]
    shared = set(rk['DEG_TSR']['drug_key']) & set(rk['DEG_sRGES']['drug_key'])
    for n, r in rk.items():
        m = r['drug_key'].isin(shared)
        r['percentile_shared_DEG'] = np.nan
        r.loc[m, 'percentile_shared_DEG'] = r.loc[m, 'score'].rank(pct=True) * 100
    print(f'drugs shared by DEG_TSR and DEG_sRGES: {len(shared)}')
    for n, r in rk.items():
        if 'n_instances' in r:
            r['stratum'] = pd.cut(r['n_instances'], STRATA, labels=STRATA_LABELS).astype(str)
            r['n_drugs_stratum'] = r.groupby('stratum')['score'].transform('size')
            r['percentile_stratum'] = r.groupby('stratum')['score'].rank(pct=True) * 100
    st = rk['DEG_sRGES'].groupby('stratum').agg(n=('score', 'size'),
                                                  median_best_pct=('drug_key', lambda k: rk['DEG_best'].set_index('drug_key')
                                                                   .loc[k, 'percentile'].median()))
    print('instance-count strata of DEG_sRGES, with the median global percentile of DEG_best:')
    print(st.reindex(STRATA_LABELS).round(1).to_string())
    print(f"Spearman(DEG_best score, n_instances) = "
          f"{spearmanr(rk['DEG_best']['score'], rk['DEG_best']['n_instances']).statistic:.3f}")

    # universe sizes
    x = md.copy()
    x['pert_dose'] = pd.to_numeric(x['pert_dose'], errors='coerce')
    x['pert_time'] = pd.to_numeric(x['pert_time'], errors='coerce')
    kept = x[(x['pert_dose'] > 0) & x['pert_time'].isin([6, 24])].assign(drug_key=lambda d: key(d['pert_iname']))
    times = kept.groupby('drug_key')['pert_time'].agg(lambda t: set(t))
    pi_keys, tsr_keys = set(rk['DEG_sRGES']['drug_key']), set(rk['DEG_TSR']['drug_key'])
    only_pi = pi_keys - tsr_keys
    only6 = {k for k in only_pi if times.get(k) == {6}}
    print(f'universe: per-instance {len(pi_keys)}, TSR {len(tsr_keys)}, shared {len(pi_keys & tsr_keys)}, '
          f'per-instance only {len(only_pi)} (6 h-only instances: {len(only6)}), TSR only {len(tsr_keys - pi_keys)}')
    print(f'TSR-only drugs: {sorted(tsr_keys - pi_keys)[:30]}')
    shared_6only = {k for k in (pi_keys & tsr_keys) if times.get(k) == {6}}
    print(f'shared drugs with 6 h-only instances in the per-instance route: {len(shared_6only)}')

    long = pd.concat([r.assign(ranking=n) for n, r in rk.items()], ignore_index=True)
    cols = ['ranking', 'drug_key', 'drug', 'score', 'rank', 'percentile', 'percentile_shared_DEG',
            'n_instances', 'n_cell_lines', 'stratum', 'n_drugs_stratum', 'percentile_stratum', 'cs_run_id']
    long = long.reindex(columns=cols)
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(P.oout_dir / f'{P.slug}_drug_rankings.tsv', sep='\t', index=False)

    kn = long[long['drug_key'].isin(KNOWN_DRUGS)].copy()
    kn['n_drugs'] = kn['ranking'].map({n: len(r) for n, r in rk.items()})
    m3 = pd.concat([r.assign(ranking=n) for n, r in mean3.items()])
    m3 = m3[m3['drug_key'].isin(KNOWN_DRUGS)].assign(n_drugs=lambda d: d['ranking'].map({n: len(r) for n, r in mean3.items()}))
    kn = pd.concat([kn, m3], ignore_index=True)
    kn = kn[['drug_key', 'ranking', 'score', 'rank', 'n_drugs', 'percentile', 'percentile_shared_DEG',
             'n_instances', 'n_cell_lines', 'stratum', 'n_drugs_stratum', 'percentile_stratum',
             'cs_run_id']].sort_values(['drug_key', 'ranking'])
    kn.to_csv(P.oout_dir / f'{P.slug}_known_drug_ranks.tsv', sep='\t', index=False)
    missing = [d for d in KNOWN_DRUGS for n in rk if d not in set(rk[n]['drug_key'])]
    print(f'known drugs missing from a ranking: {missing or "none"}')
    print(kn.round(4).to_string(index=False))

    ag = agreement(rk, list(rk))
    strata = rk['DEG_sRGES'].groupby('stratum')['drug_key'].apply(set)
    tsr_names = [n for n in rk if n.endswith('_TSR')]
    for lab in STRATA_LABELS:
        if lab in strata:
            a2 = agreement(rk, tsr_names + ['DEG_sRGES'], stratum=lab, keep=strata[lab])
            ag = pd.concat([ag, a2[a2['ranking_b'] == 'DEG_sRGES']], ignore_index=True)
    ag.to_csv(P.oout_dir / f'{P.slug}_ranking_agreement.tsv', sep='\t', index=False)
    print(ag[(ag.stratum == 'all') & (ag.ranking_a.str.startswith('DEG') | ag.ranking_b.str.startswith('DEG'))]
          .round(3).to_string(index=False))
    print(ag[(ag.stratum != 'all') & (ag.ranking_a == 'DEG_TSR')].round(3).to_string(index=False))

    inst = known_instances(md, weights)
    inst.to_csv(P.oout_dir / f'{P.slug}_known_drug_instances.tsv', sep='\t', index=False)
    print(inst.round(4).to_string(index=False))

    panel(rk, ag, kn, deg_names)
