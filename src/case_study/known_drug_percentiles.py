"""
Step 'validation', known-drug part with cc.KNOWN_DRUGS_REPORT = 'percentiles':
position of the known drugs (cc.KNOWN_DRUGS), and agreement between drug
rankings, for the drug-collapsed (TSR/LMM) connectivity-score runs and the
per-instance sRGES scores of step 'srges'.

Design
------
- Rankings (one score per drug, more negative = stronger reversal):
  * '<config>_TSR': the most recent cs_batch run of conf.DISEASE on drug
    input 'all' (drugs collapsed before CS), per configuration in
    TSR_CONFIGS, resolved with loader.load_cs_run (every configuration must
    have a run, the step stops otherwise); score = pooled 6h_24h
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
- Known drugs: cc.KNOWN_DRUGS (name: class; class labels only, no statistic
  uses them). For each, the number of L1000 instances, of instances on LINCS
  lines with a CCLE weight, and the lines without a weight that carry its
  instances.
- Lines without CCLE expression (no weight): their instances enter the
  unweighted sRGES and are left out of the weighted one; counted overall and
  per known drug, with the drugs scored by sRGES but not by wsRGES.
- The two routes use different LINCS instance pools: the '_TSR' runs come
  from GSE92742 level-3 profiles collapsed by a linear mixed model, the
  per-instance rankings from the level-5 gold set of the Chen 2017
  distribution; a comparison between them compares method and pool together.
- Instance-count strata: for the per-instance rankings the known drugs also
  get a percentile within their stratum of instance count (STRATA, number of
  instances the ranking averages: 6h / 24h with dose > 0 for sRGES_unweighted
  and best, of those the ones on lines with a weight for sRGES_weighted), and
  the Spearman correlation of DEG_TSR with each DEG per-instance ranking is
  also given per stratum.
- With cc.CS_CHECKS, the known drugs are also ranked in the per-compound
  sRGES / wsRGES recomputed from the 1-based per-instance scores of
  cs_checks ('<...>_rank1based' rankings).

Inputs
------
- P.oout_dir/<slug>_srges_per_drug.tsv, <slug>_lincs_ccle_mapping.tsv,
  <slug>_cell_line_weights.tsv
- conf.cs_log_filename and the CS files it points to
- conf.LINCS_METADATA_PATH

Outputs
-------
- P.oout_dir/<slug>_known_drug_ranks.tsv
    one row per (known drug, ranking, collapse): score, rank, n_drugs,
    percentile_own, percentile_shared, cs_run_id (TSR) and, for per-instance
    rankings, n_instances, n_cell_lines, instance_stratum,
    n_drugs_stratum, percentile_stratum.
- P.oout_dir/<slug>_known_drug_instances.tsv
    per known drug: L1000 instances (all, 6h / 24h with dose > 0, on lines
    with a weight), lines without a weight, in the TSR universe.
- P.oout_dir/<slug>_ranking_comparison.tsv
    one row per (pair of rankings, stratum): n_shared, spearman,
    top50_overlap, top100_overlap; stratum 'all' for every pair, instance-
    count strata for DEG_TSR vs the DEG per-instance rankings.
- P.img_dir/<slug>_tsr_vs_srges_panel.png
    A: percentile of the known drugs in every ranking; B: DEG_TSR vs
    DEG_sRGES percentile on the shared drugs, known drugs marked; C: Spearman
    correlation between rankings.
"""

from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

import conf
from loader import load_cs_run, read_lincs_metadata
from case_study.common import figures, colours

TSR_CONFIGS = {'DEG': (0, 0, 0), 'DEG_LM': (0, 1, 0), 'mith': (1, 0, 0), 'mith_LM': (1, 1, 0), 'mith_pw': (1, 0, 1)}
SRGES_NAMES = {'sRGES_unweighted': 'sRGES', 'sRGES_weighted': 'wsRGES', 'best': 'best'}
DEG_RANKINGS = ['DEG_TSR', 'DEG_sRGES', 'DEG_wsRGES', 'DEG_best']
TIMEPOINTS = ('6h_24h', '24h', '6h')
TOP_K = (50, 100)
STRATA = ([0, 1, 2, 5, 9, 25, 50, np.inf], ['1', '2', '3-5', '6-9', '10-25', '26-50', '>50'])


def key(s):
    return pd.Series(s, dtype=str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True).values


def load_tsr():
    """{'<config>_TSR': DataFrame(drug_key, pooled, best_timepoint)}, {label: cs_run_id}."""
    out, ids = {}, {}
    for c, (m, lm, pw) in TSR_CONFIGS.items():
        try:
            df, path, run_id = load_cs_run(disease=conf.DISEASE, landmark_disease=False, drug_run_id='all', mith=m,
                                           cs_on_LM=lm, CS_ON_PATHWAYS=pw, CS_METHOD=conf.CS_METHOD, verbose=False)
        except FileNotFoundError as e:
            raise FileNotFoundError(f"known-drug report 'percentiles' needs a drug-collapsed (drug input 'all') "
                                    f"CS run of {conf.DISEASE} for every configuration; none for {c}_TSR "
                                    f"(mith={m}, cs_on_LM={lm}, CS_ON_PATHWAYS={pw})") from e
        tp = df['pert_id'].str.extract(rf'^(.*?)_({"|".join(TIMEPOINTS)})$')
        df['drug_key'] = key(tp[0])
        df['tp'] = tp[1]
        pooled = df[df.tp == '6h_24h'].groupby('drug_key')['connectivity_score'].mean()
        best = df.groupby('drug_key')['connectivity_score'].min()
        out[f'{c}_TSR'] = pd.DataFrame({'pooled': pooled, 'best_timepoint': best})
        ids[f'{c}_TSR'] = run_id
        print(f'{c}_TSR: {run_id}; {df.drug_key.nunique()} drugs, {len(df)} rows')
    return out, ids


def load_per_instance(P):
    t = pd.read_csv(P.oout_dir / f'{P.slug}_srges_per_drug.tsv', sep='\t')
    t['drug_key'] = key(t['pert_iname'])
    out = {}
    for (c, col), g in t.groupby(['config', 'collapse']):
        g = g.groupby('drug_key').agg(score=('score', 'mean'), n_instances=('n_instances', 'sum'),
                                       n_cell_lines=('n_cell_lines', 'max'))
        out[f'{c}_{SRGES_NAMES[col]}'] = g
    return out


def percentile(s):
    return s.rank(method='average') / s.notna().sum() * 100


def panel(cc, known, pv, order, R, comp):
    P = cc.P
    BLUE, RED = colours(cc)
    fig = plt.figure(figsize=(16, 11))
    gsp = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], wspace=0.35)
    a = fig.add_subplot(gsp[0, :])
    cols = plt.cm.tab10.colors
    for j, name in enumerate(order):
        sub = pv[pv.ranking == name].set_index('drug')
        y = np.arange(len(known)) + (j - len(order) / 2) * 0.045
        x = [sub.loc[d, 'percentile_own'] if d in sub.index else np.nan for d in known]
        mk = 's' if name.endswith('_TSR') else ('o' if name.endswith('sRGES') else '^')
        a.scatter(x, y, marker=mk, s=30, color=cols[j % 10], label=name,
                  edgecolor='k' if name.startswith('DEG') else 'none', linewidth=.5)
    a.set_yticks(np.arange(len(known)))
    a.set_yticklabels([f'{d}\n({c.split(",")[0]})' for d, c in known.items()], fontsize=9)
    a.set_xlim(0, 100)
    a.invert_yaxis()
    a.set_xlabel('percentile within the ranking (low = stronger reversal)')
    a.set_title(f'A. Known {figures(cc)["known_drugs_label"]} drugs: percentile in each ranking '
                '(squares _TSR, circles sRGES/wsRGES, triangles best; missing = not scored)')
    a.legend(fontsize=8, ncol=2, loc='center left', bbox_to_anchor=(1.0, .5))

    b = fig.add_subplot(gsp[1, 0])
    if 'DEG_sRGES' in R:
        s = pd.concat([R['DEG_TSR'], R['DEG_sRGES']], axis=1, join='inner', keys=['t', 's'])
        pt, ps = percentile(s.t), percentile(s.s)
        b.scatter(pt, ps, s=4, color='0.6', alpha=.5)
        for d in known:
            k = key([d])[0]
            if k in s.index:
                b.scatter(pt[k], ps[k], s=40, color=RED)
                b.annotate(d, (pt[k], ps[k]), xytext=(4, 4), textcoords='offset points', fontsize=9)
        rho = spearmanr(s.t, s.s).statistic
        b.set_xlabel('DEG_TSR percentile (pooled 6h_24h CS)')
        b.set_ylabel('DEG_sRGES percentile (unweighted sRGES)')
        b.set_title(f'B. DEG_TSR vs DEG_sRGES on {len(s)} shared drugs (Spearman {rho:.2f})')

    c = fig.add_subplot(gsp[1, 1])
    M = pd.DataFrame(np.eye(len(order)), index=order, columns=order)
    for _, r in comp[comp.stratum == 'all'].iterrows():
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


def run(cc, rank1based=None):
    """rank1based: the 1-based per-compound rankings returned by cs_checks.run, or None."""
    P = cc.P
    known = cc.KNOWN_DRUGS
    md = read_lincs_metadata(dtype={'id': str})
    mapping = pd.read_csv(P.oout_dir / f'{P.slug}_lincs_ccle_mapping.tsv', sep='\t')
    no_w = set(mapping.loc[~mapping.has_expression, 'lincs_cell_id'])

    tsr, tsr_ids = load_tsr()
    inst = load_per_instance(P)
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
    for d, cls_ in known.items():
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

    if rank1based:
        extra = []
        for name, s in rank1based.items():
            for d, cls_ in known.items():
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

    panel(cc, known, pv, order, R, comp)
