"""
Connectivity-score Monte Carlo p-value filter of a case study: drug rankings
before and after dropping the connectivity scores (CS) that are not
significant against random gene sets.

Design
------
- p-value of a CS: column cs_p_value of a CS run, the fraction of N_ITER = 1000
  random gene sets of the same size whose |CS| is >= the observed |CS|
  (no pseudocount, so p = 0 occurs). It is two-sided in the sign of the CS: a
  significant CS can be negative (reversal) or positive.
- Filter: p' = (k + 1) / (N_ITER + 1) with k = cs_p_value x N_ITER;
  p_adj = Benjamini-Hochberg adjustment of p' over the rows of one CS run;
  'bh05' keeps p_adj <= ALPHA (primary), 'raw05' keeps p' <= ALPHA
  (secondary). With n0 rows at k = 0 in a run of m rows, the p_adj of those
  rows is at most (1 / (N_ITER + 1)) x m / n0 (the BH step-up takes the
  minimum over later ranks, so it can be lower).
- Offset variants: 'bh05' rebuilds the dose/time offset table from the surviving
  instances; 'bh05f' (frozen) scores the surviving instances with the offset
  table of the unfiltered instances (sRGES.build_reference_diff is replaced
  by that table while per_drug_scores.collapse runs).
- Per-instance route (one CS per LINCS instance; 71 cell lines x configurations
  of step 'srges'): the filter is applied to every CS run (one run = one
  configuration x one cell line, all its instances) and the instances that do
  not pass are dropped; per_drug_scores.collapse then recomputes sRGES, wsRGES,
  'best' and the top-line variants from the surviving instances (the dose/time
  correction table is rebuilt from them); a compound with no surviving
  instance leaves the ranking.
- TSR route (one CS per drug, drug input 'all', pooled 6h_24h rows): p_adj is
  computed over the 6h_24h rows of the run; drugs that do not pass leave the
  ranking.
- Rankings, ranks, percentiles, shared-universe percentiles and instance-count
  strata are built with the functions of known_drugs, with the same columns
  as <slug>_drug_rankings.tsv.
- Sign convention: more negative = stronger reversal of the disease signature.

Inputs
------
- <slug>_srges_cs_runs_used.tsv, <slug>_cell_line_weights.tsv,
  <slug>_top_cell_lines.tsv in P.oout_dir
- the CS runs resolved with loader.load_cs_run (per-instance: by run ID;
  TSR: drug input 'all', configurations of known_drugs.TSR_CONFIGS)
- conf.LINCS_METADATA_PATH, conf.srges_diff_mode

Usage
-----
    import sys; sys.path.insert(0, '<repo>/src')
    from case_study import cs_pvalue_filter as pf
    res = pf.compute(cc)      # cc = case_study_pipeline.load_case_conf(slug)
"""

from contextlib import contextmanager

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import conf
import sRGES
from loader import load_cs_run, read_lincs_metadata
from case_study import per_drug_scores as pds
from case_study import known_drugs as kd

N_ITER = 1000
ALPHA = 0.05
TOP = 20
MODES = ('none', 'bh05', 'raw05', 'bh05f')
LONG_COLS = ['ranking', 'config', 'collapse', 'drug_key', 'drug', 'score', 'rank', 'percentile',
             'percentile_shared_DEG', 'n_instances', 'n_cell_lines', 'stratum', 'n_drugs_stratum',
             'percentile_stratum', 'cs_run_id']


# ------------------------------------------------------------ statistics

def bh(p):
    """Benjamini-Hochberg adjusted p-values (same order as p)."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p, kind='mergesort')
    adj = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(adj, 1.0)
    return out


def flag(p):
    """p' = (k + 1) / (N_ITER + 1) of cs_p_value p, its BH adjustment, and the two keep masks."""
    k = np.rint(np.asarray(p, dtype=float) * N_ITER)
    pp = (k + 1) / (N_ITER + 1)
    padj = bh(pp)
    return pd.DataFrame({'p_prime': pp, 'p_adj': padj, 'keep_bh05': padj <= ALPHA, 'keep_raw05': pp <= ALPHA})


def run_stats(df):
    """Resolution of one CS run: rows, rows at p = 0, p_adj of the rows at p = 0 (exact, and its upper bound
    (1 / (N_ITER + 1)) x m / n0), kept rows."""
    f = flag(df['cs_p_value'])
    z = df['cs_p_value'].values == 0
    n0 = int(z.sum())
    m = len(df)
    return {'n_rows': m, 'n_p0': n0, 'frac_p0': n0 / m if m else np.nan,
            'p_adj_at_p0': float(f.loc[z, 'p_adj'].max()) if n0 else np.nan,
            'p_adj_p0_upper_bound': (1 / (N_ITER + 1)) * m / n0 if n0 else np.nan,
            'n_keep_bh05': int(f['keep_bh05'].sum()), 'n_keep_raw05': int(f['keep_raw05'].sum())}


def md_table(df, digits=3):
    """Markdown pipe table of a DataFrame: floats rounded to `digits` (scientific notation below
    10^-digits), columns named n_*, rank* or top* printed as integers."""
    def fmt(v, c):
        if isinstance(v, (float, np.floating)):
            if np.isnan(v):
                return ''
            if str(c).startswith(('n_', 'rank', 'top')) and float(v).is_integer():
                return str(int(v))
            return f'{v:.2e}' if v != 0 and abs(v) < 10 ** -digits else f'{v:.{digits}f}'
        return str(v)
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '|'.join(['---'] * len(cols)) + '|']
    for _, r in df.iterrows():
        lines.append('| ' + ' | '.join(fmt(r[c], c) for c in cols) + ' |')
    return '\n'.join(lines)


# ----------------------------------------------------------------- inputs

def load_inputs(cc):
    """Instance metadata, cell-line weights, top-line sets and the CS runs used by step 'srges'."""
    P = cc.P
    md = read_lincs_metadata(dtype={'id': str})[['id', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose', 'pert_time']]
    w = pd.read_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t')
    weights = dict(zip(w['lincs_cell_id'], w['cor']))
    t = pd.read_csv(P.oout_dir / f'{P.slug}_top_cell_lines.tsv', sep='\t')
    top = {cut: g['lincs_cell_id'].dropna().tolist() for cut, g in t.groupby('cutoff', sort=False)}
    top = {cut: lines for cut, lines in top.items() if lines}
    used = pd.read_csv(P.oout_dir / f'{P.slug}_srges_cs_runs_used.tsv', sep='\t')
    return md, weights, top, used


def load_instances(slug, used, config):
    """All instances of one configuration with their CS, p-values and keep masks (flag per CS run)."""
    parts, stats = [], []
    for _, r in used[(used['config'] == config) & used['cs_run_id'].notna()].iterrows():
        df, _, run_id = load_cs_run(disease=slug, landmark_disease=False, selected_cs_run_id=r['cs_run_id'],
                                    verbose=False)
        df = df[['LINCS_id', 'connectivity_score', 'cs_p_value']].reset_index(drop=True)
        f = flag(df['cs_p_value'])
        df = pd.concat([df, f], axis=1).assign(cell_line=r['cell_line'], cs_run_id=run_id)
        parts.append(df)
        stats.append({'route': 'per_instance', 'config': config, 'cell_line': r['cell_line'],
                      'cs_run_id': run_id, **run_stats(df)})
    return pd.concat(parts, ignore_index=True), pd.DataFrame(stats)


def load_tsr(slug):
    """Pooled 6h_24h rows of the drug-collapsed CS runs, per '<config>_TSR', with p-values and keep masks."""
    rows, stats = {}, []
    for label, (m, lm, pw) in kd.TSR_CONFIGS.items():
        try:
            df, _, run_id = load_cs_run(disease=slug, landmark_disease=False, drug_run_id='all', mith=m, cs_on_LM=lm,
                                        CS_ON_PATHWAYS=pw, CS_METHOD='bin_chen', verbose=False)
        except FileNotFoundError:
            continue
        dt = df['pert_id'].astype(str).map(kd.split_pert_id)
        df = df.assign(drug=dt.str[0], tp=dt.str[1])
        df = df[df['tp'] == '6h_24h'].reset_index(drop=True)
        df['drug_key'] = kd.key(df['drug'])
        df = pd.concat([df[['drug', 'drug_key', 'connectivity_score', 'cs_p_value']], flag(df['cs_p_value'])], axis=1)
        df['cs_run_id'] = run_id
        rows[f'{label}_TSR'] = df
        stats.append({'route': 'tsr', 'config': label, 'cell_line': 'all drugs', 'cs_run_id': run_id,
                      **run_stats(df)})
    return rows, pd.DataFrame(stats)


# --------------------------------------------------------------- rankings

def tsr_ranking(rows, keep=None):
    """Drug ranking of one TSR run (rows of the surviving drugs when `keep` is a mask), as known_drugs does."""
    r = rows if keep is None else rows[keep]
    g = r.groupby('drug_key').agg(score=('connectivity_score', 'mean'), drug=('drug', 'first')).reset_index()
    g['cs_run_id'] = rows['cs_run_id'].iloc[0]
    return g


@contextmanager
def frozen_offsets(diff):
    """Make sRGES.compute_sRGES use the dose/time offset table `diff` instead of building one."""
    orig = sRGES.build_reference_diff
    sRGES.build_reference_diff = lambda *a, **k: diff
    try:
        yield
    finally:
        sRGES.build_reference_diff = orig


def instance_table(inst, md, weights, top, keep=None, frozen=None):
    """Per-drug scores of every collapse (per_drug_scores.collapse) from all instances or the kept ones;
    `frozen` = offset table (per_drug_scores.offsets of the unfiltered instances) to use instead of a rebuilt one."""
    assert conf.srges_diff_mode == 'corrected', conf.srges_diff_mode
    x = (inst if keep is None else inst[keep])[['LINCS_id', 'connectivity_score']]
    if frozen is None:
        return pds.collapse(x, md, weights, top)
    with frozen_offsets(frozen):
        return pds.collapse(x, md, weights, top)


def line_stats(inst, md, weights, config, mode):
    """Sign split of the kept CS and cell-line retention of one configuration under filter `mode`.

    Instances used by sRGES = dose > 0 at 6 or 24 h. Weight mass = sum of the cell-line weights of the used
    instances on weighted lines (instance level), or of the weighted lines having a used instance (line level)."""
    x = inst.assign(LINCS_id=inst['LINCS_id'].astype(str), keep=inst[f'keep_{mode}'].values)
    x = x.merge(md[['id', 'cell_id', 'pert_dose', 'pert_time']], left_on='LINCS_id', right_on='id')
    u = x[(pd.to_numeric(x['pert_dose'], errors='coerce') > 0) &
          pd.to_numeric(x['pert_time'], errors='coerce').isin([6, 24])].copy()
    u['w'] = u['cell_id'].map(weights)
    uw = u.dropna(subset=['w'])
    lw_all = uw.groupby('cell_id')['w'].first()
    lw_kept = uw[uw['keep']].groupby('cell_id')['w'].first()
    cs = x.loc[x['keep'], 'connectivity_score']
    return {'config': config, 'mode': mode, 'n_kept': int(x['keep'].sum()), 'n_kept_neg': int((cs < 0).sum()),
            'n_kept_pos': int((cs > 0).sum()), 'n_used': len(u), 'n_used_kept': int(u['keep'].sum()),
            'n_lines_used': int(u['cell_id'].nunique()), 'n_lines_kept': int(u.loc[u['keep'], 'cell_id'].nunique()),
            'n_weighted_lines': len(lw_all), 'n_weighted_lines_kept': len(lw_kept),
            'frac_line_weight_lost': 1 - lw_kept.sum() / lw_all.sum(),
            'frac_instance_weight_lost': 1 - uw.loc[uw['keep'], 'w'].sum() / uw['w'].sum()}


def build_rankings(tsr, pi):
    """Rankings and their long table from TSR rankings {name: DataFrame} and per-instance tables
    (DataFrame with config, collapse, pert_iname, score, n_instances, n_cell_lines)."""
    s = pi.assign(drug_key=kd.key(pi['pert_iname']))
    rk, meta = {}, {}
    for n, r in tsr.items():
        rk[n] = kd.add_ranks(r)
        meta[n] = (n[:-len('_TSR')], 'pooled')
    pi_names = []
    for (config, coll), d in s.groupby(['config', 'collapse']):
        n = f'{config}_{kd.SCORE_TYPES[coll]}'
        r = d.groupby('drug_key').agg(score=('score', 'mean'), drug=('pert_iname', 'first'),
                                      n_instances=('n_instances', 'sum'),
                                      n_cell_lines=('n_cell_lines', 'max')).reset_index()
        rk[n] = kd.add_ranks(r)
        meta[n] = (config, kd.SCORE_TYPES[coll])
        pi_names.append(n)
    shared = (set(rk['DEG_TSR']['drug_key']) & set(rk['DEG_sRGES']['drug_key'])) \
        if 'DEG_TSR' in rk and 'DEG_sRGES' in rk else set()
    for r in rk.values():
        kd.add_shared_percentile(r, shared)
    for n in pi_names:
        kd.add_strata(rk[n])
    cfgs = list(kd.TSR_CONFIGS) + sorted({meta[n][0] for n in pi_names} - set(kd.TSR_CONFIGS))
    order = [n for c in cfgs for n in [f'{c}_TSR'] + [f'{c}_{t}' for t in kd.SCORE_TYPES.values()] if n in rk]
    long = pd.concat([rk[n].assign(ranking=n, config=meta[n][0], collapse=meta[n][1]) for n in order],
                     ignore_index=True).reindex(columns=LONG_COLS)
    return {n: rk[n] for n in order}, long, meta


def compute(cc):
    """Unfiltered ('none'), BH-filtered ('bh05'), raw-filtered ('raw05') and BH-filtered with frozen offsets
    ('bh05f') rankings of a case study.

    Returns a dict: 'rk' {mode: {ranking: DataFrame}}, 'long' {mode: DataFrame}, 'meta' {ranking: (config, collapse)},
    'run_stats' (one row per CS run), 'pvals' {(route, config): array of cs_p_value}, 'inst_kept' (per
    configuration: instances before / after, per mode)."""
    slug = cc.P.slug
    md, weights, top, used = load_inputs(cc)
    tsr_rows, st_tsr = load_tsr(slug)
    stats = [st_tsr]
    frozen_check, offs = [], {}
    pvals = {('tsr', n[:-4]): r['cs_p_value'].values for n, r in tsr_rows.items()}
    tsr_by_mode = {m: {n: tsr_ranking(r, None if m == 'none' else r[f'keep_{m.rstrip("f")}']) for n, r in tsr_rows.items()}
                   for m in MODES}
    tabs = {m: [] for m in MODES}
    inst_kept, lstats = [], []
    for config in pds.CONFIGS:
        inst, st = load_instances(slug, used, config)
        stats.append(st)
        pvals[('per_instance', config)] = inst['cs_p_value'].values
        frozen = pds.offsets(inst[['LINCS_id', 'connectivity_score']], md)
        for m in MODES:
            keep = None if m == 'none' else inst[f'keep_{m.rstrip("f")}']
            t = instance_table(inst, md, weights, top, keep, frozen if m == 'bh05f' else None)
            tabs[m].append(t.assign(config=config))
            inst_kept.append({'config': config, 'mode': m, 'n_instances': len(inst),
                              'n_kept': len(inst) if keep is None else int(keep.sum())})
        for m in ('bh05', 'raw05'):
            lstats.append(line_stats(inst, md, weights, config, m))
        # frozen offsets with every instance reproduce the unfiltered scores
        t0 = instance_table(inst, md, weights, top, None, frozen)
        t1 = tabs['none'][-1].drop(columns='config')
        j = t0.merge(t1, on=['collapse', 'pert_iname'], suffixes=('_f', '_u'))
        frozen_check.append({'config': config, 'n_compared': len(j), 'max_abs_diff': float((j['score_f'] - j['score_u']).abs().max())})
        offs[config] = frozen.to_dict()
    out = {'rk': {}, 'long': {}, 'tables': {}}
    for m in MODES:
        pi = pd.concat(tabs[m], ignore_index=True)
        out['tables'][m] = pi
        out['rk'][m], out['long'][m], meta_m = build_rankings(tsr_by_mode[m], pi)
        if m == 'none':
            meta = meta_m   # every ranking of the unfiltered build
    out.update({'meta': meta, 'run_stats': pd.concat(stats, ignore_index=True), 'pvals': pvals,
                'inst_kept': pd.DataFrame(inst_kept), 'tsr_rows': tsr_rows, 'line_stats': pd.DataFrame(lstats),
                'frozen_check': pd.DataFrame(frozen_check), 'offsets': offs})
    return out


# ------------------------------------------------------------ comparisons

def sanity(long_new, stored, tol=1e-12):
    """Compare a recomputed long ranking table with a stored one.

    Per column: rows compared, largest absolute difference, rows differing by more than `tol`, and of those
    the rows that are not near-ties. A near-tie is a drug whose stored score lies within `tol` of another
    drug's score in the same ranking: ranks and percentiles of such drugs depend on floating-point rounding
    (about 1e-16) of the two scores, so they may differ by one or two positions. A column passes when every
    difference above `tol` is on a near-tie. Also returns the number of (ranking, drug) rows in only one table."""
    a = long_new.set_index(['ranking', 'drug_key'])
    b = stored.set_index(['ranking', 'drug_key'])
    only_a, only_b = len(a.index.difference(b.index)), len(b.index.difference(a.index))
    j = a.join(b, how='inner', lsuffix='_new', rsuffix='_old')
    near = pd.Series(False, index=j.index)
    for n, g in j.groupby(level='ranking'):
        d = g['score_old'].sort_values().diff().abs()
        f = (d <= tol) | (d.shift(-1) <= tol)
        near.loc[f.index[f.values]] = True
    rows = []
    for c in ['score', 'rank', 'percentile', 'percentile_shared_DEG', 'n_instances', 'n_cell_lines',
              'n_drugs_stratum', 'percentile_stratum']:
        x, y = j[f'{c}_new'].astype(float), j[f'{c}_old'].astype(float)
        both = x.notna() & y.notna()
        diff = (x - y).abs()
        bad = both & (diff > tol)
        rows.append({'column': c, 'n_compared': int(both.sum()),
                     'max_abs_diff': float(diff[both].max()) if both.any() else np.nan,
                     'n_diff_above_tol': int(bad.sum()), 'n_diff_not_near_tie': int((bad & ~near).sum()),
                     'nan_mismatch': int((x.isna() != y.isna()).sum())})
    r = pd.DataFrame(rows)
    r['pass'] = (r['n_diff_not_near_tie'] == 0) & (r['nan_mismatch'] == 0)
    r.loc[r['column'] == 'score', 'pass'] &= r['max_abs_diff'] <= tol
    return r, only_a, only_b


def compare_rankings(rk_a, rk_b, k=TOP):
    """Per ranking: universe sizes, Spearman correlation of the scores on the shared drugs and
    top-k overlap (own universes, and restricted to the shared drugs)."""
    rows = []
    for n, a in rk_a.items():
        b = rk_b.get(n)
        row = {'ranking': n, 'n_before': len(a), 'n_after': 0 if b is None else len(b)}
        if b is not None and len(b) > 2:
            x, y = a.set_index('drug_key')['score'], b.set_index('drug_key')['score']
            sh = x.index.intersection(y.index)
            row.update({'n_shared': len(sh),
                        'spearman': spearmanr(x[sh], y[sh]).statistic if len(sh) > 2 else np.nan,
                        f'top{k}_overlap': len(set(x.nsmallest(k).index) & set(y.nsmallest(k).index)),
                        f'top{k}_overlap_shared': len(set(x[sh].nsmallest(k).index) & set(y[sh].nsmallest(k).index))})
        rows.append(row)
    r = pd.DataFrame(rows)
    r['frac_kept'] = r['n_after'] / r['n_before']
    return r


def known_table(rk_a, rk_b, known, mode):
    """Position of the known drugs before (rk_a) and after (rk_b) the filter.

    known: {display name: (class, [candidate drug keys])}; the first candidate present in a ranking is used."""
    def find(r, keys):
        if r is None:
            return None
        for k in keys:
            h = r[r['drug_key'] == k]
            if len(h):
                return h.iloc[0]
        return None
    rows = []
    for n, a in rk_a.items():
        b = rk_b.get(n)
        for d, (cls, keys) in known.items():
            ha, hb = find(a, keys), find(b, keys)
            rows.append({'filter': mode, 'ranking': n, 'drug': d, 'class': cls,
                         'n_before': len(a), 'rank_before': np.nan if ha is None else ha['rank'],
                         'percentile_before': np.nan if ha is None else ha['percentile'],
                         'n_after': 0 if b is None else len(b),
                         'rank_after': np.nan if hb is None else hb['rank'],
                         'percentile_after': np.nan if hb is None else hb['percentile'],
                         'status': 'not scored' if ha is None else ('kept' if hb is not None else 'removed')})
    t = pd.DataFrame(rows)
    t['delta_percentile'] = t['percentile_after'] - t['percentile_before']
    return t
