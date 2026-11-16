"""
Step 'validation', known-drug report (run when cc.KNOWN_DRUGS is declared):
drug rankings from the drug-collapsed (TSR/LMM) connectivity-score runs and
from the per-instance sRGES scores of step 'srges', position of the known
drugs (cc.KNOWN_DRUGS, name: class) in every ranking, agreement between
rankings, per-instance CS of the known drugs.

Design
------
- Rankings (one score per drug, more negative = stronger reversal):
  * '<config>_TSR': the pooled 6h_24h connectivity_score of the most recent
    cs_batch run of conf.DISEASE on drug input 'all' (drugs collapsed by the
    TSR/LMM model before CS), per configuration in TSR_CONFIGS, resolved with
    loader.load_cs_run. A configuration without such a run is skipped with
    a message; DEG_TSR is required.
  * '<config>_sRGES', '<config>_wsRGES', '<config>_best': the
    'sRGES_unweighted', 'sRGES_weighted' and 'best' per-drug scores of step
    'srges', for every configuration present in its output; DEG_sRGES is
    required.
  * '<config>_sRGES_topQ3', '_wsRGES_topQ3', '_sRGES_topMed',
    '_wsRGES_topMed': the top-line sRGES collapses of step 'srges'
    ('sRGES_(un)weighted_top_q3' / '_top_median'), when present.
  * Two more collapses of every TSR run, reported for the known drugs only:
    '<config>_TSR_mean3tp' = mean of the drug's 6h, 24h and 6h_24h scores,
    '<config>_TSR_besttp' = the most negative of them.
  * With cc.CS_CHECKS, the per-compound sRGES / wsRGES recomputed by
    cs_checks from 1-based per-instance scores
    ('<config>_sRGES_rank1based', '<config>_wsRGES_rank1based'), reported
    for the known drugs only.
- Drug key: the L1000 compound name in lower case, alphanumerics only
  (pert_iname in the per-instance route, the drug part of pert_id in the TSR
  route); names that collapse to the same key take the mean score (per-
  instance: summed n_instances, largest n_cell_lines).
- Rank: ascending score, ties at the lowest rank. Percentile:
  rank(pct=True) x 100 of the score (average rank of ties / n; low
  percentile = strong reverser):
  * 'percentile': over the full universe of that ranking;
  * 'percentile_shared_DEG': over the drugs scored by both DEG_TSR and
    DEG_sRGES (for a ranking that scores only part of them, e.g.
    DEG_wsRGES, over that part);
  * 'percentile_stratum' (per-instance rankings): among the drugs of the
    same ranking and instance-count stratum. The stratum bins n_instances,
    the number of instances that ranking averages (6 h / 24 h with dose > 0;
    for wsRGES, of those the ones on lines with a weight), by STRATA.
- Instance pools: the TSR route collapses the LINCS GSE92742 level-3
  instances of a drug (per-drug LMM of tsr/); the per-instance route uses the
  level-5 gold set of conf.LINCS_METADATA_PATH. A TSR vs per-instance
  comparison therefore differs in both the collapse and the instances.
- Agreement of two rankings on their shared drugs: Spearman correlation and
  overlap of the top-k (k in TOP_K), the top-k taken after restricting both
  rankings to the shared drugs; for every pair of rankings (stratum 'all'),
  and for every TSR ranking vs every per-instance ranking also within each
  instance-count stratum of the per-instance ranking.
- Per-instance CS of the known drugs: every L1000 instance of a known drug in
  conf.LINCS_METADATA_PATH, with its CS in the run of its cell line listed in
  <slug>_srges_cs_runs_used.tsv, for every configuration scored by step
  'srges'; 'used_by_srges' = dose > 0 and time 6 or 24 h (the instances
  sRGES and best keep); 'weight' = cell-line weight of step 'weights' (empty
  = line without CCLE expression, left out of wsRGES).

Inputs
------
- P.oout_dir/<slug>_srges_per_drug.tsv, P.oout_dir/<slug>_srges_cs_runs_used.tsv
- P.oout_dir/<slug>_cell_line_weights.tsv, P.oout_dir/<slug>_lincs_ccle_mapping.tsv
- conf.cs_log_filename and the CS files it points to
- conf.LINCS_METADATA_PATH

Outputs
-------
- P.oout_dir/<slug>_drug_rankings.tsv
    one row per (ranking, drug_key) of the main rankings: config, collapse,
    drug, score, rank, percentile, percentile_shared_DEG, n_instances,
    n_cell_lines, stratum, n_drugs_stratum, percentile_stratum (per-instance
    only), cs_run_id (TSR only).
- P.oout_dir/<slug>_known_drug_ranks.tsv
    one row per (known drug, ranking), main rankings plus the extra TSR
    collapses and the 1-based rankings: drug, drug_key, class, config,
    collapse, score, rank, n_drugs, percentile, percentile_shared_DEG,
    n_instances, n_cell_lines, stratum, n_drugs_stratum, percentile_stratum,
    cs_run_id; empty score = drug not scored by that ranking.
- P.oout_dir/<slug>_ranking_agreement.tsv
    one row per (pair of rankings, stratum): ranking_a, ranking_b, stratum
    ('all' = no stratification), n_a, n_b (drugs of each ranking), n_shared,
    spearman, top<k>_overlap for every k in TOP_K.
- P.oout_dir/<slug>_known_drug_instances.tsv
    one row per (config, L1000 instance of a known drug): drug, class,
    pert_iname, cell_id, pert_id, id, pert_dose, pert_time, used_by_srges,
    weight, connectivity_score.
- P.img_dir/<slug>_tsr_vs_srges_panel.png
    A: percentile of the known drugs in every ranking (colour =
    configuration, marker = score type); B: DEG_TSR vs DEG_sRGES percentile
    on the shared drugs, known drugs marked; C: Spearman correlation between
    rankings (all shared drugs); D: per-instance CS of the known drugs per
    configuration (instances used by sRGES).
"""

from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

import conf
from loader import load_cs_run, read_lincs_metadata
from case_study.common import figures, colours, config_colour

TSR_CONFIGS = {'DEG': (0, 0, 0), 'DEG_LM': (0, 1, 0), 'mith': (1, 0, 0), 'mith_LM': (1, 1, 0), 'mith_pw': (1, 0, 1)}
SCORE_TYPES = {'sRGES_unweighted': 'sRGES', 'sRGES_weighted': 'wsRGES', 'best': 'best',
               'sRGES_unweighted_top_q3': 'sRGES_topQ3', 'sRGES_weighted_top_q3': 'wsRGES_topQ3',
               'sRGES_unweighted_top_median': 'sRGES_topMed', 'sRGES_weighted_top_median': 'wsRGES_topMed'}
MARKERS = {'TSR': 's', 'sRGES': 'o', 'wsRGES': 'D', 'best': '^',
           'sRGES_topQ3': 'v', 'wsRGES_topQ3': 'P', 'sRGES_topMed': '<', 'wsRGES_topMed': 'X'}
TOP_K = [50, 100]
STRATA = [0, 1, 2, 5, 9, 25, 50, np.inf]
STRATA_LABELS = ['1', '2', '3-5', '6-9', '10-25', '26-50', '>50']
TIMEPOINTS = ('6h_24h', '24h', '6h')


def key(s):
    return pd.Series(s, dtype=str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True).values


def split_pert_id(pid):
    for tp in TIMEPOINTS:
        if pid.endswith('_' + tp):
            return pid[:-(len(tp) + 1)], tp
    return pid, None


# -------------------------------------------------------------- rankings

def tsr_rankings():
    """Main '<config>_TSR' rankings (pooled 6h_24h) and the extra collapses, with (config, collapse) of each."""
    main, extra, meta = {}, {}, {}
    for label, (m, lm, pw) in TSR_CONFIGS.items():
        try:
            df, path, run_id = load_cs_run(disease=conf.DISEASE, landmark_disease=False, drug_run_id='all', mith=m,
                                           cs_on_LM=lm, CS_ON_PATHWAYS=pw, CS_METHOD=conf.CS_METHOD, verbose=False)
        except FileNotFoundError as e:
            if label == 'DEG':
                raise FileNotFoundError(f"known-drug report: no drug-collapsed (drug input 'all') CS run of "
                                        f"{conf.DISEASE} for DEG_TSR (mith=0, cs_on_LM=0, CS_ON_PATHWAYS=0); "
                                        f"DEG_TSR is required") from e
            print(f"{label}_TSR: no drug-collapsed (drug input 'all') CS run of {conf.DISEASE} "
                  f"(mith={m}, cs_on_LM={lm}, CS_ON_PATHWAYS={pw}); configuration skipped")
            continue
        dt = df['pert_id'].astype(str).map(split_pert_id)
        df = df.assign(drug=dt.str[0], tp=dt.str[1])
        df['drug_key'] = key(df['drug'])
        g = df[df['tp'] == '6h_24h'].groupby('drug_key')
        main[f'{label}_TSR'] = g.agg(score=('connectivity_score', 'mean'), drug=('drug', 'first')).reset_index()
        meta[f'{label}_TSR'] = (label, 'pooled')
        g = df.groupby('drug_key')
        extra[f'{label}_TSR_mean3tp'] = g.agg(score=('connectivity_score', 'mean'), drug=('drug', 'first')).reset_index()
        extra[f'{label}_TSR_besttp'] = g.agg(score=('connectivity_score', 'min'), drug=('drug', 'first')).reset_index()
        meta[f'{label}_TSR_mean3tp'] = (label, 'mean3tp')
        meta[f'{label}_TSR_besttp'] = (label, 'best_timepoint')
        for n in (f'{label}_TSR', f'{label}_TSR_mean3tp', f'{label}_TSR_besttp'):
            (main if n in main else extra)[n]['cs_run_id'] = run_id
        print(f'{label}_TSR: {run_id}; {len(main[f"{label}_TSR"])} drugs (6h_24h), '
              f'{len(extra[f"{label}_TSR_mean3tp"])} drugs (any timepoint)')
    return main, extra, meta


def per_instance_rankings(P):
    s = pd.read_csv(P.oout_dir / f'{P.slug}_srges_per_drug.tsv', sep='\t')
    s['drug_key'] = key(s['pert_iname'])
    out, meta = {}, {}
    for (config, coll), d in s.groupby(['config', 'collapse']):
        r = d.groupby('drug_key').agg(score=('score', 'mean'), drug=('pert_iname', 'first'),
                                      n_instances=('n_instances', 'sum'),
                                      n_cell_lines=('n_cell_lines', 'max')).reset_index()
        out[f'{config}_{SCORE_TYPES[coll]}'] = r
        meta[f'{config}_{SCORE_TYPES[coll]}'] = (config, SCORE_TYPES[coll])
    return out, meta


def add_ranks(r):
    r = r.copy()
    r['rank'] = r['score'].rank(method='min')
    r['percentile'] = r['score'].rank(pct=True) * 100
    return r


def add_shared_percentile(r, shared):
    m = r['drug_key'].isin(shared)
    r['percentile_shared_DEG'] = np.nan
    r.loc[m, 'percentile_shared_DEG'] = r.loc[m, 'score'].rank(pct=True) * 100


def add_strata(r):
    r['stratum'] = pd.cut(r['n_instances'], STRATA, labels=STRATA_LABELS).astype(str)
    r['n_drugs_stratum'] = r.groupby('stratum')['score'].transform('size')
    r['percentile_stratum'] = r.groupby('stratum')['score'].rank(pct=True) * 100


def agreement(rk, pairs, stratum='all', keep=None):
    rows = []
    for a, b in pairs:
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


# ------------------------------------------------------ known instances

def known_instances(P, md, weights, known):
    """Every L1000 instance of a known drug with its CS in every scored configuration."""
    used = pd.read_csv(P.oout_dir / f'{P.slug}_srges_cs_runs_used.tsv', sep='\t')
    configs = [c for c, u in used.groupby('config', sort=False) if u['cs_run_id'].notna().all()]
    ids = md[md['drug_key'].isin(known['drug_key'])].copy()
    ids['drug'] = ids['drug_key'].map(dict(zip(known['drug_key'], known.index)))
    ids['class'] = ids['drug_key'].map(dict(zip(known['drug_key'], known['class'])))
    dose = pd.to_numeric(ids['pert_dose'], errors='coerce')
    time = pd.to_numeric(ids['pert_time'], errors='coerce')
    ids['used_by_srges'] = (dose > 0) & time.isin([6, 24])
    ids['weight'] = ids['cell_id'].map(weights)
    rows = []
    for config in configs:
        runs = used[used['config'] == config].set_index('cell_line')['cs_run_id']
        for cell, sub in ids.groupby('cell_id'):
            df, _, _ = load_cs_run(selected_cs_run_id=runs[cell], disease=conf.DISEASE, landmark_disease=False,
                                   verbose=False)
            cs = df.assign(LINCS_id=df['LINCS_id'].astype(str)).set_index('LINCS_id')['connectivity_score']
            rows.append(sub.assign(config=config, connectivity_score=sub['id'].map(cs)))
    x = pd.concat(rows, ignore_index=True)
    n_miss = x['connectivity_score'].isna().sum()
    if n_miss:
        print(f'known-drug instances without a CS in the run of their line: {n_miss}')
    cols = ['config', 'drug', 'class', 'pert_iname', 'cell_id', 'pert_id', 'id', 'pert_dose', 'pert_time',
            'used_by_srges', 'weight', 'connectivity_score']
    x['config'] = pd.Categorical(x['config'], configs)
    x['drug'] = pd.Categorical(x['drug'], list(known.index))
    return x[cols].sort_values(['config', 'drug', 'connectivity_score']).reset_index(drop=True)


# ----------------------------------------------------------------- panel

def panel(cc, rk, meta, order, ag, kn, inst, known):
    P = cc.P
    BLUE, RED = colours()
    F = figures(cc)
    nk = len(known)
    fig = plt.figure(figsize=(18, 13 + max(4.5, 0.9 * nk + 2)))
    gsp = fig.add_gridspec(3, 2, height_ratios=[max(4.5, 0.9 * nk + 2), 7.5, 5.5], width_ratios=[1, 1.15],
                           wspace=0.35, hspace=0.45)

    # A: percentile of the known drugs in every ranking
    a = fig.add_subplot(gsp[0, :])
    step = 0.8 / max(len(order), 1)
    for j, name in enumerate(order):
        config, typ = meta[name]
        typ = 'TSR' if typ == 'pooled' else typ
        sub = kn[kn['ranking'] == name].set_index('drug')['percentile']
        y = np.arange(nk) - 0.4 + step * (j + 0.5)
        a.scatter(sub.reindex(known.index).values, y, marker=MARKERS[typ], s=34,
                  color=config_colour(config), edgecolor='k', linewidth=.4, zorder=3)
    a.axvline(50, color='k', lw=.6)
    a.set_yticks(np.arange(nk))
    a.set_yticklabels([f'{d}\n({c})' for d, c in known['class'].items()], fontsize=9)
    a.set_ylim(nk - 0.5, -0.5)
    a.set_xlim(0, 100)
    a.set_xlabel('percentile within the ranking (low = stronger reversal)')
    a.set_title(f'A. Known {F["known_drugs_label"]} drugs: percentile in every ranking (own universe; '
                f'missing = not scored)')
    configs = list(dict.fromkeys(meta[n][0] for n in order))
    h1 = [Line2D([], [], marker='o', ls='', color=config_colour(c), markeredgecolor='k', label=c) for c in configs]
    h2 = [Line2D([], [], marker=m, ls='', color='0.6', markeredgecolor='k', label=t) for t, m in MARKERS.items()
          if t in {'TSR' if meta[n][1] == 'pooled' else meta[n][1] for n in order}]
    l1 = a.legend(handles=h1, title='configuration', loc='upper left', bbox_to_anchor=(1.0, 1.0))
    a.add_artist(l1)
    a.legend(handles=h2, title='score', loc='lower left', bbox_to_anchor=(1.0, 0.0))

    # B: DEG_TSR vs DEG_sRGES on the shared drugs
    b = fig.add_subplot(gsp[1, 0])
    x = rk['DEG_TSR'].set_index('drug_key')['percentile_shared_DEG'].dropna()
    y = rk['DEG_sRGES'].set_index('drug_key')['percentile_shared_DEG'].dropna()
    sh = x.index.intersection(y.index)
    b.scatter(x[sh], y[sh], s=4, color='0.6', alpha=.5, rasterized=True)
    for j, (d, k) in enumerate(known['drug_key'].items()):
        if k in sh:
            b.scatter(x[k], y[k], s=45, color=RED, edgecolor='k', zorder=3)
            b.annotate(d, (x[k], y[k]), xytext=(6, 6 - 14 * (j % 3)), textcoords='offset points', fontsize=9,
                       arrowprops=dict(arrowstyle='-', color='0.3', lw=.5))
    b.plot([0, 100], [0, 100], color='k', lw=.6)
    rho = spearmanr(x[sh], y[sh]).statistic
    b.set_xlim(0, 100); b.set_ylim(0, 100)
    b.set_xlabel('DEG_TSR percentile (pooled 6h_24h CS)')
    b.set_ylabel('DEG_sRGES percentile (unweighted sRGES)')
    b.set_title(f'B. DEG_TSR vs DEG_sRGES on {len(sh):,} shared drugs (Spearman {rho:.2f})')

    # C: Spearman between rankings
    c = fig.add_subplot(gsp[1, 1])
    M = pd.DataFrame(np.eye(len(order)), index=order, columns=order)
    for _, r in ag[ag['stratum'] == 'all'].iterrows():
        M.loc[r.ranking_a, r.ranking_b] = M.loc[r.ranking_b, r.ranking_a] = r.spearman
    im = c.imshow(M.values, vmin=-1, vmax=1, cmap='bwr')
    c.set_xticks(range(len(order)))
    c.set_xticklabels(order, rotation=60, ha='right', fontsize=8)
    c.set_yticks(range(len(order)))
    c.set_yticklabels(order, fontsize=8)
    c.tick_params(axis='y', pad=2)
    for i in range(len(order)):
        for j in range(len(order)):
            c.text(j, i, f'{M.values[i, j]:.1f}', ha='center', va='center', fontsize=8)
    c.grid(False)
    c.set_title('C. Spearman correlation between rankings (shared drugs)')
    fig.colorbar(im, ax=c, fraction=.046)

    # D: per-instance CS of the known drugs
    d = fig.add_subplot(gsp[2, :])
    u = inst[inst['used_by_srges']]
    cfgs = list(u['config'].cat.categories) if hasattr(u['config'], 'cat') else list(dict.fromkeys(u['config']))
    w = 0.8 / max(len(cfgs), 1)
    rng = np.random.default_rng(0)
    for i, drug in enumerate(known.index):
        for j, cfg in enumerate(cfgs):
            v = u.loc[(u['drug'] == drug) & (u['config'] == cfg), 'connectivity_score'].values
            x0 = i - 0.4 + w * (j + 0.5)
            d.scatter(x0 + rng.uniform(-w / 4, w / 4, len(v)), v, s=14, color=config_colour(cfg),
                      edgecolor='k', linewidth=.3, zorder=3)
            if len(v):
                d.plot([x0 - w / 2.5, x0 + w / 2.5], [np.median(v)] * 2, color='k', lw=1.2, zorder=4)
    d.axhline(0, color='k', lw=.6)
    d.set_xticks(np.arange(nk))
    d.set_xticklabels(known.index)
    d.set_xlim(-0.5, nk - 0.5)
    d.set_ylabel('per-instance connectivity score')
    d.set_title('D. Per-instance CS of the known drugs (instances with dose > 0, 6 or 24 h; bar = median)')
    d.legend(handles=[Line2D([], [], marker='o', ls='', color=config_colour(c), markeredgecolor='k', label=c)
                      for c in cfgs], title='configuration', loc='upper left', bbox_to_anchor=(1.0, 1.0))

    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_tsr_vs_srges_panel.png')
    plt.close(fig)


# ------------------------------------------------------------------- run

def run(cc, rank1based=None):
    """rank1based: the 1-based per-compound rankings returned by cs_checks.run, or None."""
    P = cc.P
    known = pd.DataFrame({'class': pd.Series(cc.KNOWN_DRUGS)})
    known.index.name = 'drug'
    known['drug_key'] = key(known.index)
    md = read_lincs_metadata(dtype={'id': str})
    md = md[['id', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose', 'pert_time']].copy()
    md['drug_key'] = key(md['pert_iname'])
    w = pd.read_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t')
    weights = dict(zip(w['lincs_cell_id'], w['cor']))
    mapping = pd.read_csv(P.oout_dir / f'{P.slug}_lincs_ccle_mapping.tsv', sep='\t')
    no_w = set(mapping.loc[~mapping['has_expression'], 'lincs_cell_id'])

    tsr, tsr_extra, meta = tsr_rankings()
    pi, meta_pi = per_instance_rankings(P)
    meta.update(meta_pi)
    if 'DEG_sRGES' not in pi:
        raise ValueError(f'known-drug report: {P.slug}_srges_per_drug.tsv has no DEG sRGES_unweighted scores '
                         f'(DEG_sRGES is required)')
    rk = {n: add_ranks(r) for n, r in {**tsr, **pi}.items()}
    extra = {n: add_ranks(r) for n, r in tsr_extra.items()}
    order = [n for c in list(TSR_CONFIGS) + sorted({m[0] for m in meta_pi.values()} - set(TSR_CONFIGS))
             for n in [f'{c}_TSR'] + [f'{c}_{t}' for t in SCORE_TYPES.values()] if n in rk]
    print('rankings: ' + ', '.join(f'{n} ({len(rk[n])})' for n in order))

    shared = set(rk['DEG_TSR']['drug_key']) & set(rk['DEG_sRGES']['drug_key'])
    for r in list(rk.values()) + list(extra.values()):
        add_shared_percentile(r, shared)
    print(f'drugs shared by DEG_TSR and DEG_sRGES: {len(shared)}')
    for n in pi:
        add_strata(rk[n])

    # strata and universes
    st = rk['DEG_sRGES'].groupby('stratum').agg(
        n=('score', 'size'),
        median_best_pct=('drug_key', lambda k: rk['DEG_best'].set_index('drug_key').loc[k, 'percentile'].median()))
    print('instance-count strata of DEG_sRGES, with the median global percentile of DEG_best:')
    print(st.reindex(STRATA_LABELS).round(1).to_string())
    print(f"Spearman(DEG_best score, n_instances) = "
          f"{spearmanr(rk['DEG_best']['score'], rk['DEG_best']['n_instances']).statistic:.3f}")
    un = set(rk['DEG_TSR']['drug_key']) - set(md['drug_key'])
    print(f'DEG_TSR drugs whose name is not an L1000 pert_iname of the metadata: {len(un)}'
          + (f' (e.g. {", ".join(sorted(un)[:8])})' if un else ''))
    x = md.copy()
    x['pert_dose'] = pd.to_numeric(x['pert_dose'], errors='coerce')
    x['pert_time'] = pd.to_numeric(x['pert_time'], errors='coerce')
    kept = x[(x['pert_dose'] > 0) & x['pert_time'].isin([6, 24])]
    times = kept.groupby('drug_key')['pert_time'].agg(lambda t: set(t))
    pi_keys, tsr_keys = set(rk['DEG_sRGES']['drug_key']), set(rk['DEG_TSR']['drug_key'])
    only_pi = pi_keys - tsr_keys
    only6 = {k for k in only_pi if times.get(k) == {6}}
    print(f'universe: per-instance {len(pi_keys)}, TSR {len(tsr_keys)}, shared {len(pi_keys & tsr_keys)}, '
          f'per-instance only {len(only_pi)} (6 h-only instances: {len(only6)}), TSR only {len(tsr_keys - pi_keys)}')
    print(f'TSR-only drugs: {sorted(tsr_keys - pi_keys)[:30]}')
    print(f'shared drugs with 6 h-only instances in the per-instance route: '
          f'{len({k for k in (pi_keys & tsr_keys) if times.get(k) == {6}})}')
    n_nw = kept['cell_id'].isin(no_w).sum()
    print(f'6h/24h instances with dose > 0: {len(kept)}; on the {len(no_w)} lines without a weight: {n_nw} '
          f'({n_nw / len(kept):.1%})')
    if 'DEG_wsRGES' in rk:
        lost = pi_keys - set(rk['DEG_wsRGES']['drug_key'])
        print(f'drugs scored by DEG_sRGES but not by DEG_wsRGES (every instance on a line without weight): {len(lost)}')

    # drug rankings
    long = pd.concat([r.assign(ranking=n, config=meta[n][0], collapse=meta[n][1]) for n, r in
                      ((n, rk[n]) for n in order)], ignore_index=True)
    cols = ['ranking', 'config', 'collapse', 'drug_key', 'drug', 'score', 'rank', 'percentile',
            'percentile_shared_DEG', 'n_instances', 'n_cell_lines', 'stratum', 'n_drugs_stratum',
            'percentile_stratum', 'cs_run_id']
    long = long.reindex(columns=cols)
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(P.oout_dir / f'{P.slug}_drug_rankings.tsv', sep='\t', index=False)

    # known drugs in every ranking
    allr = {n: rk[n] for n in order}
    allr.update(extra)
    for n, s in (rank1based or {}).items():
        allr[n] = add_ranks(pd.DataFrame({'drug_key': s.index, 'score': s.values}))
        meta[n] = (n.split('_sRGES')[0].split('_wsRGES')[0],
                   'wsRGES_rank1based' if '_wsRGES_' in n else 'sRGES_rank1based')
    kcols = ['score', 'rank', 'percentile', 'percentile_shared_DEG', 'n_instances', 'n_cell_lines', 'stratum',
             'n_drugs_stratum', 'percentile_stratum', 'cs_run_id']
    krows = []
    for d, kk in known['drug_key'].items():
        for n, r in allr.items():
            rec = {'drug': d, 'drug_key': kk, 'class': known.loc[d, 'class'], 'ranking': n,
                   'config': meta[n][0], 'collapse': meta[n][1], 'n_drugs': len(r)}
            hit = r[r['drug_key'] == kk]
            if len(hit):
                rec.update({c: hit.iloc[0][c] for c in kcols if c in r})
            krows.append(rec)
    kn = pd.DataFrame(krows).reindex(columns=['drug', 'drug_key', 'class', 'ranking', 'config', 'collapse',
                                              'score', 'rank', 'n_drugs'] + kcols[2:])
    kn.to_csv(P.oout_dir / f'{P.slug}_known_drug_ranks.tsv', sep='\t', index=False)
    missing = sorted({f'{d} ({n})' for d, n in kn.loc[kn['score'].isna() & kn['ranking'].isin(order),
                                                       ['drug', 'ranking']].values})
    print(f'known drugs missing from a ranking: {missing or "none"}')
    piv = kn[kn['ranking'].isin(order)].pivot_table(index='drug', columns='ranking', values='percentile',
                                                    dropna=False)
    print('\nknown drugs, percentile (own universe):')
    print(piv.reindex(index=known.index, columns=order).round(1).to_string())
    print('\nknown drugs, percentile_shared_DEG:')
    print(kn[kn['ranking'].isin(order)].pivot_table(index='drug', columns='ranking', values='percentile_shared_DEG',
                                                    dropna=False).reindex(index=known.index, columns=order)
          .round(1).to_string())
    print('\nknown drugs, percentile within the instance-count stratum (per-instance rankings):')
    print(kn[kn['percentile_stratum'].notna()].pivot_table(index=['drug', 'stratum'], columns='ranking',
                                                           values='percentile_stratum').round(1).to_string())
    print('\nknown drugs, percentile in the extra TSR collapses and 1-based rankings:')
    print(kn[~kn['ranking'].isin(order)].pivot_table(index='drug', columns='ranking', values='percentile')
          .round(1).to_string())

    # agreement
    pairs = list(combinations(order, 2))
    ag = [agreement(rk, pairs)]
    tsr_names = [n for n in order if meta[n][1] == 'pooled']
    pi_names = [n for n in order if n in pi]
    for b in pi_names:
        strata = rk[b].groupby('stratum')['drug_key'].apply(set)
        for lab in STRATA_LABELS:
            if lab in strata:
                ag.append(agreement(rk, [(t, b) for t in tsr_names], stratum=lab, keep=strata[lab]))
    ag = pd.concat(ag, ignore_index=True)
    ag.to_csv(P.oout_dir / f'{P.slug}_ranking_agreement.tsv', sep='\t', index=False)
    print('\nranking agreement, DEG_TSR vs the other rankings (shared drugs):')
    print(ag[(ag['stratum'] == 'all') & (ag['ranking_a'] == 'DEG_TSR')].round(3).to_string(index=False))
    print('\n<config>_TSR vs the per-instance rankings of the same configuration, by instance-count stratum (Spearman):')
    same = ag[(ag['stratum'] != 'all') & (ag['ranking_a'].map(lambda n: meta[n][0]) ==
                                          ag['ranking_b'].map(lambda n: meta[n][0]))]
    print(same.pivot_table(index='stratum', columns='ranking_b', values='spearman')
          .reindex(STRATA_LABELS).round(3).to_string())

    # per-instance CS of the known drugs
    inst = known_instances(P, md, weights, known)
    inst.to_csv(P.oout_dir / f'{P.slug}_known_drug_instances.tsv', sep='\t', index=False)
    first = inst[inst['config'] == inst['config'].cat.categories[0]]
    print('\nknown drugs, L1000 instances:')
    for d in known.index:
        s = first[first['drug'] == d]
        k = s[s['used_by_srges']]
        print(f'  {d}: {len(s)} instances, {len(k)} with dose > 0 at 6 / 24 h '
              f'({int(k["weight"].notna().sum())} on lines with a weight); times {sorted(int(t) for t in k["pert_time"].unique())}; '
              f'lines without weight {sorted(set(k["cell_id"]) & no_w)}; in the DEG_TSR universe: '
              f'{known.loc[d, "drug_key"] in tsr_keys}')
    print(inst[inst['used_by_srges']].groupby(['config', 'drug'], observed=True)['connectivity_score']
          .agg(['size', 'median', 'min']).round(4).to_string())

    panel(cc, rk, meta, order, ag, kn, inst, known)
