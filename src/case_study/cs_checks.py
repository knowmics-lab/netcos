"""
Step 'validation', CS checks (run when cc.CS_CHECKS is declared):
known-answer checks of the per-instance DEG connectivity scores of the study.

Design
------
(i) For CS_CHECKS['n_check'] random instances of CS_CHECKS['cell_line']
    (seed 0), the stored connectivity_score is recomputed with
    connectivity_score.bin_chen_connectivity on the same inputs, and with an
    independent RGES on 1-based positions.
(ii) For every instance of every cell line, an independent implementation of
    the Chen 2017 RGES (cmap_score_new: KS statistics a, b of the up- and
    down-regulated disease genes on the positions of the drug ranking,
    ks = a if a > b else -b, score ks_up - ks_down, also when the two have
    the same sign) is computed with the 1-based positions of the published
    code and with the 0-based positions of connectivity_score.rank_genes
    (a + 1/r, b - 1/r), both on the drug ranking ordered as in rank_genes
    (descending pandas sort, so tied drug values keep the same order). The
    0-based score is compared with the stored CS; the 1-based score is
    collapsed per compound with sRGES.compute_sRGES (unweighted and weighted)
    to measure the effect of the position convention on the per-compound
    ranking. The share of instances whose ks_up and ks_down have the same
    sign is reported (Lamb 2006 scores those instances 0; Chen 2017 and this
    pipeline do not).

Stored CS: the most recent DEG cs_batch run of conf.DISEASE on drug input
'<cell_id>_LM' (loader.load_cs_run).

Inputs
------
- CS_CHECKS['disease_signature']: the DEG disease signature of the CS runs
  (TSR format: ';' separator, ',' decimal)
- the per-instance DEG drug inputs (conf.CS_IN_DRUG_DIR/<cell_id>_LM)
- conf.cs_log_filename and the CS files it points to, conf.LINCS_METADATA_PATH
- P.oout_dir/<slug>_cell_line_weights.tsv

Outputs
-------
- P.oout_dir/<slug>_cs_known_answer_check.tsv
    per checked instance of (i): stored CS, CS recomputed with the pipeline
    function, CS of the independent implementation, ks_up, ks_down.
- P.oout_dir/<slug>_cs_position_convention.tsv
    per instance of (ii): LINCS_id, cell_id, cs_stored, cs_0based,
    cs_1based, ks_up_0based, ks_down_0based.
- returned to the known-drug report: per-compound DEG sRGES / wsRGES from the
  1-based per-instance scores ('DEG_sRGES_rank1based', 'DEG_wsRGES_rank1based').
"""

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import spearmanr

import conf
from loader import load_cs_run, load_single_signature_cs_input, standardize_disease_signature, read_lincs_metadata
from connectivity_score import get_common_genes, bin_chen_connectivity
from sRGES import compute_sRGES


def key(s):
    return pd.Series(s, dtype=str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True).values


def drug_dir(cell):
    return conf.CS_IN_DRUG_DIR / (cell + '_LM')


def read_signature(C):
    dis = pd.read_csv(C['disease_signature'], sep=';', decimal=',', dtype={'gene_id': str})
    return standardize_disease_signature(dis, source='tsr')


def disease_common(C):
    """DEG disease signature restricted to the genes of the drug signatures."""
    dis = read_signature(C)
    first = next(drug_dir(C['cell_line']).glob('*_signature_gene_id.pkl'))
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


def position_convention_check(C, md, weights):
    dis = disease_common(C)
    up = dis.loc[dis.DE_log2_FC > 0, 'gene_id'].values
    down = dis.loc[dis.DE_log2_FC < 0, 'gene_id'].values
    jobs = []
    for c in sorted(md.cell_id.unique()):
        df, _, _ = load_cs_run(disease=conf.DISEASE, landmark_disease=False, drug_run_id=c + '_LM', mith=0, cs_on_LM=0,
                               CS_ON_PATHWAYS=0, CS_METHOD=conf.CS_METHOD, verbose=False)
        jobs.append((c, list(zip(df['LINCS_id'].astype(str), df['connectivity_score']))))
    res = Parallel(n_jobs=C['n_jobs'])(delayed(ks_one_line)(c, ids, up, down) for c, ids in jobs)
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


def known_answer_check(C, md):
    cell = C['cell_line']
    m = md[md.cell_id == cell]
    df, path, run_id = load_cs_run(disease=conf.DISEASE, landmark_disease=False, drug_run_id=cell + '_LM',
                                   mith=0, cs_on_LM=0, CS_ON_PATHWAYS=0, CS_METHOD=conf.CS_METHOD, verbose=False)
    df['LINCS_id'] = df['LINCS_id'].astype(str)
    ids = df['LINCS_id'].sample(n=min(C['n_check'], len(df)), random_state=0).tolist()
    dis = read_signature(C)
    ddir = drug_dir(cell)
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
    print(f'\nknown-answer check, {cell}_LM ({run_id}), {len(ch)} instances, '
          f'{len(dis_c)} disease genes on the drug genes ({(dis_c.DE_log2_FC > 0).sum()} up, '
          f'{(dis_c.DE_log2_FC < 0).sum()} down), r = {r}:')
    print(f'  max |stored - pipeline recomputed| = {d1:.2e}')
    print(f'  max |stored - independent RGES|    = {d2:.2e}  (1/r = {1 / r:.2e}); '
          f'Spearman {spearmanr(ch.cs_stored, ch.cs_independent).statistic:.4f}')
    print(f'  instances with ks_up and ks_down of the same sign: {same:.0%}')
    return ch


def run(cc):
    """Write both check tables; return the 1-based per-compound rankings of (ii)."""
    P = cc.P
    C = cc.CS_CHECKS
    md = read_lincs_metadata(dtype={'id': str})
    w = pd.read_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t')
    weights = dict(zip(w['lincs_cell_id'], w['cor']))
    ch = known_answer_check(C, md)
    conv, r1 = position_convention_check(C, md, weights)
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    ch.to_csv(P.oout_dir / f'{P.slug}_cs_known_answer_check.tsv', sep='\t', index=False)
    conv.to_csv(P.oout_dir / f'{P.slug}_cs_position_convention.tsv', sep='\t', index=False)
    return r1
