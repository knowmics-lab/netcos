"""
Step 'validation', CS checks (run when cc.CS_CHECKS is declared):
recomputation checks of the per-instance connectivity scores of the study,
for every configuration scored by step 'srges' (DEG, mith, mith_LM,
mith_pw).

Design
------
- Stored CS: for every configuration and LINCS cell line, the cs_batch run
  listed in P.oout_dir/<slug>_srges_cs_runs_used.tsv (the runs step 'srges'
  used), opened with loader.load_cs_run. A configuration is checked only
  when every cell line has a run there.
- Inputs of a run: its row of conf.cs_log_filename (disease_run_id,
  drug_input_dir, disease_input_dir, mith, cs_on_LM, CS_ON_PATHWAYS). The
  disease and drug signatures are read as cs_batch reads them
  (loader.load_disease_signature, loader.load_single_signature_cs_input):
  DEG = DE_log2_FC per gene, mith / mith_LM = MITHrIL Perturbation per node,
  mith_pw = MITHrIL Corrected Accumulator per pathway; with cs_on_LM both
  are restricted to the LINCS landmark genes; the disease signature is then
  restricted to the items of the drug signature.
(i) Recomputation check: for n_check random instances of CS_CHECKS['cell_line']
    (seed 0) and every configuration, the stored connectivity_score is
    recomputed with cs_batch.run_connectivity_score_drugs_batch (the
    function of the run) on the run's inputs, and with an independent RGES
    (see (ii)).
(ii) Independent RGES on every instance: for every instance of every cell
    line and every configuration, an independent implementation of the
    Chen 2017 RGES (cmap_score_new: KS statistics a, b of the up- and
    down-regulated disease items on their positions 1..r in the drug
    ranking, ks = a if a > b else -b, score ks_up - ks_down, also when the
    two have the same sign), on the drug ranking ordered as in rank_genes
    (descending pandas sort, so tied drug values keep the same order), is
    compared with the stored CS (overall and per cell line), leaving out
    the instances where a KS statistic is a tie (|a - b| < TIE for ks_up or
    ks_down: the choice ks = a or -b is then decided by floating-point
    rounding and differs between equivalent implementations; their number is
    reported and they are flagged in 'ks_tie'). The share of instances
    whose ks_up and ks_down have the same sign is reported (Lamb 2006
    scores those instances 0; Chen 2017 and this pipeline do not).
- CS_CHECKS = {'cell_line': ..., 'n_check': ..., 'n_jobs': ...}; n_check
  (instances of (i), default N_CHECK) and n_jobs (parallel jobs of (ii),
  default N_JOBS) are optional.

Inputs
------
- P.oout_dir/<slug>_srges_cs_runs_used.tsv
- conf.cs_log_filename and the CS files and input folders it points to
- conf.LINCS_METADATA_PATH, conf.LINCS_BC_DATA (landmark genes)

Outputs
-------
- P.oout_dir/<slug>_cs_recomputation_check.tsv
    per (config, checked instance) of (i): cs_run_id, LINCS_id, pert_iname,
    cs_stored, cs_pipeline_recomputed, cs_independent, ks_up, ks_down.
- P.oout_dir/<slug>_cs_independent_rges.tsv
    per (config, instance) of (ii): LINCS_id, cell_id, cs_stored,
    cs_independent, ks_up, ks_down, ks_tie.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import spearmanr

import conf
from loader import (load_cs_run, load_disease_signature, load_single_signature_cs_input,
                    load_landmark_gene_ids, read_lincs_metadata)
from connectivity_score import get_common_genes
from cs_batch import run_connectivity_score_drugs_batch
from case_study.per_drug_scores import CONFIGS

N_CHECK = 20   # instances recomputed in (i)
N_JOBS = 12    # parallel jobs of (ii)
CHUNK = 500    # instances per parallel job of (ii)
TIE = 1e-12    # |a - b| below which a KS statistic is a tie (ks = a or -b decided by rounding)


def key(s):
    return pd.Series(s, dtype=str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True).values


# ------------------------------------------------------------ run inputs

def columns(mith, pw):
    """(item id column, value column, p-value column) of a configuration, as in cs_batch."""
    if pw:
        return 'ID', 'Corrected Accumulator', 'Adjusted pValue'
    return 'gene_id', ('Perturbation' if mith else 'DE_log2_FC'), 'adj.p.value'


def run_record(run_id):
    """Row of conf.cs_log_filename of a CS run of conf.DISEASE."""
    log = pd.read_csv(conf.cs_log_filename, sep='\t', dtype=str)
    r = log[(log['cs_run_id'] == run_id) & (log['disease_run_id'] == conf.DISEASE)]
    if r.empty:
        raise ValueError(f'{run_id}: no row of disease {conf.DISEASE} in {conf.cs_log_filename}')
    r = r.iloc[-1]
    return {'run_id': run_id, 'disease_run_id': r['disease_run_id'], 'drug_dir': r['drug_input_dir'],
            'disease_dir': r['disease_input_dir'], 'mith': int(r['mith']), 'lm': int(r['cs_on_LM']),
            'pw': int(float(r['CS_ON_PATHWAYS']))}


def load_drug(i, rec, lm_ids):
    """Drug signature of instance i as cs_batch uses it (item id, value, p-value columns)."""
    idc, val, pv = columns(rec['mith'], rec['pw'])
    d = load_single_signature_cs_input(i, rec['drug_dir'], mith=rec['mith'], pathway=rec['pw'])
    if rec['lm']:
        d = d[d[idc].isin(lm_ids)].reset_index(drop=True)
    return d[[idc, val, pv]]


def disease_common(rec, lm_ids, first_drug):
    """Disease signature of the run restricted to the items of the drug signature."""
    idc, val, pv = columns(rec['mith'], rec['pw'])
    dis = load_disease_signature(rec['disease_run_id'], cs_input_dir=Path(rec['disease_dir']),
                                 mith=rec['mith'], pathway=rec['pw'])[[idc, val, pv]]
    if rec['lm']:
        dis = dis[dis[idc].isin(lm_ids)].reset_index(drop=True)
    i_dis, _ = get_common_genes(dis, first_drug, id_col=idc)
    return dis.iloc[i_dis].reset_index(drop=True)


# ----------------------------------------------------- independent RGES

def _ab(v, r):
    """KS statistics a, b of the sorted positions v (1..r) among r items."""
    s = len(v)
    j = np.arange(1, s + 1)
    return np.max(j / s - v / r), np.max(v / r - (j - 1) / s)


def _ks(a, b):
    return a if a > b else -b


def positions(drug, idc, val):
    """Positions 1..r of the drug items (1 = most up-regulated), ordered as connectivity_score.rank_genes orders them."""
    order = drug.sort_values(by=val, ascending=False)[idc].values
    return pd.Series(np.arange(1, len(order) + 1), index=order)


def ks_one_line(cell, ids_cs, rec, up, down, lm_ids):
    idc, val, _ = columns(rec['mith'], rec['pw'])
    rows = []
    for i, cs in ids_cs:
        d = load_drug(i, rec, lm_ids)
        r = len(d)
        pos = positions(d, idc, val)
        (au, bu), (ad, bd) = (_ab(np.sort(pos.loc[g].values), r) for g in (up, down))
        ku, kd = _ks(au, bu), _ks(ad, bd)
        rows.append({'LINCS_id': i, 'cell_id': cell, 'cs_stored': cs, 'cs_independent': ku - kd,
                     'ks_up': ku, 'ks_down': kd, 'ks_tie': abs(au - bu) < TIE or abs(ad - bd) < TIE})
    return rows


# ---------------------------------------------------------------- checks

def recomputation_check(label, run_id, C, md, lm_ids):
    cell = C['cell_line']
    rec = run_record(run_id)
    idc, val, _ = columns(rec['mith'], rec['pw'])
    df, _, _ = load_cs_run(selected_cs_run_id=run_id, disease=conf.DISEASE, landmark_disease=False, verbose=False)
    df['LINCS_id'] = df['LINCS_id'].astype(str)
    ids = df['LINCS_id'].sample(n=min(C['n_check'], len(df)), random_state=0).tolist()
    pipe, _ = run_connectivity_score_drugs_batch(
        rec['disease_run_id'], rec['mith'], ids, 0, 0, len(ids), conf.CS_METHOD, rec['drug_dir'],
        Path(rec['disease_dir']), conf.LINCS_BC_DATA, None, 1,
        save_file=False, cs_on_LM=bool(rec['lm']), cs_on_pathways=bool(rec['pw']), drug_collapsed=False)
    pipe = pipe.assign(LINCS_id=pipe['LINCS_id'].astype(str)).set_index('LINCS_id')['connectivity_score']
    dis_c = disease_common(rec, lm_ids, load_drug(ids[0], rec, lm_ids))
    names = md.set_index('id')['pert_iname']
    stored = df.set_index('LINCS_id')['connectivity_score']
    rows = []
    for i in ids:
        drug = load_drug(i, rec, lm_ids)
        # independent RGES: positions 1..r of the disease items in the drug ranking (most up = 1)
        pos = pd.Series(drug[val].rank(ascending=False, method='first').values, index=drug[idc].values)
        r = len(drug)
        ks = {}
        for side, sel in (('up', dis_c[val] > 0), ('down', dis_c[val] < 0)):
            v = np.sort(pos.loc[dis_c.loc[sel, idc]].values)
            ks[side] = _ks(*_ab(v, r))
        rows.append({'config': label, 'cs_run_id': run_id, 'LINCS_id': i, 'pert_iname': names.get(i),
                     'cs_stored': stored[i], 'cs_pipeline_recomputed': pipe[i],
                     'cs_independent': ks['up'] - ks['down'], 'ks_up': ks['up'], 'ks_down': ks['down']})
    ch = pd.DataFrame(rows)
    d1 = (ch.cs_stored - ch.cs_pipeline_recomputed).abs().max()
    d2 = (ch.cs_stored - ch.cs_independent).abs().max()
    same = (np.sign(ch.ks_up) == np.sign(ch.ks_down)).mean()
    print(f'\n{label}: recomputation check, {cell} ({run_id}), {len(ch)} instances, '
          f'{len(dis_c)} disease items on the drug items ({(dis_c[val] > 0).sum()} up, '
          f'{(dis_c[val] < 0).sum()} down), r = {r}:')
    print(f'  max |stored - pipeline recomputed| = {d1:.2e}')
    print(f'  max |stored - independent RGES|    = {d2:.2e}  (1/r = {1 / r:.2e}); '
          f'Spearman {spearmanr(ch.cs_stored, ch.cs_independent).statistic:.4f}')
    print(f'  instances with ks_up and ks_down of the same sign: {same:.0%}')
    return ch


def independent_rges_check(label, runs, C, lm_ids):
    """runs: {cell line: cs_run_id} of the configuration."""
    jobs, recs = [], []
    for c, run_id in sorted(runs.items()):
        rec = run_record(run_id)
        df, _, _ = load_cs_run(selected_cs_run_id=run_id, disease=conf.DISEASE, landmark_disease=False, verbose=False)
        jobs.append((c, rec, list(zip(df['LINCS_id'].astype(str), df['connectivity_score']))))
        recs.append(rec)
    if len({(r['disease_run_id'], r['disease_dir'], r['mith'], r['lm'], r['pw']) for r in recs}) != 1:
        raise ValueError(f'{label}: the per-instance runs differ in disease input or configuration')
    rec0 = recs[0]
    idc, val, _ = columns(rec0['mith'], rec0['pw'])
    c0, r0, ids0 = jobs[0]
    dis = disease_common(rec0, lm_ids, load_drug(ids0[0][0], r0, lm_ids))
    up = dis.loc[dis[val] > 0, idc].values
    down = dis.loc[dis[val] < 0, idc].values
    chunks = [(c, rec, ids[k:k + CHUNK]) for c, rec, ids in jobs for k in range(0, len(ids), CHUNK)]
    res = Parallel(n_jobs=C['n_jobs'])(delayed(ks_one_line)(c, ids, rec, up, down, lm_ids) for c, rec, ids in chunks)
    ch = pd.DataFrame([r for rr in res for r in rr])
    d = (ch.cs_stored - ch.cs_independent).abs().where(~ch.ks_tie, 0.0)
    print(f'\n{label}: independent RGES, all {len(ch)} instances ({len(up)} up, {len(down)} down disease items):')
    print(f'  KS ties (|a - b| < {TIE:g} for ks_up or ks_down, branch decided by rounding; excluded below): '
          f'{int(ch.ks_tie.sum())}')
    print(f'  non-tie instances: max |stored - independent| = {d.max():.2e}; with |diff| > 1e-12: {(d > 1e-12).sum()}')
    bad = ch[d > 1e-12].assign(diff=d[d > 1e-12]).groupby('cell_id')['diff'].agg(['size', 'max'])
    if len(bad):
        print('  cell lines with |stored - independent| > 1e-12 (instances, max |diff|):')
        print(bad.to_string())
    print(f'  instances with ks_up and ks_down of the same sign: '
          f'{(np.sign(ch.ks_up) == np.sign(ch.ks_down)).mean():.1%}')
    ch.insert(0, 'config', label)
    return ch


def run(cc):
    """Write both check tables."""
    P = cc.P
    C = {'n_check': N_CHECK, 'n_jobs': N_JOBS, **cc.CS_CHECKS}
    md = read_lincs_metadata(dtype={'id': str})
    md = md[['id', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose', 'pert_time']]
    used = pd.read_csv(P.oout_dir / f'{P.slug}_srges_cs_runs_used.tsv', sep='\t')
    lm_ids = [str(x) for x in load_landmark_gene_ids(conf.LINCS_BC_DATA)]
    if C['cell_line'] not in set(md['cell_id']):
        raise ValueError(f"CS_CHECKS['cell_line'] = {C['cell_line']!r} is not a LINCS cell line")

    known, conv = [], []
    for label in CONFIGS:
        u = used[used['config'] == label]
        if u.empty or u['cs_run_id'].isna().any():
            print(f'\n{label}: per-instance runs incomplete in {P.slug}_srges_cs_runs_used.tsv, not checked')
            continue
        runs = dict(zip(u['cell_line'], u['cs_run_id']))
        known.append(recomputation_check(label, runs[C['cell_line']], C, md, lm_ids))
        conv.append(independent_rges_check(label, runs, C, lm_ids))
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    if known:
        pd.concat(known, ignore_index=True).to_csv(P.oout_dir / f'{P.slug}_cs_recomputation_check.tsv',
                                                   sep='\t', index=False)
        pd.concat(conv, ignore_index=True).to_csv(P.oout_dir / f'{P.slug}_cs_independent_rges.tsv',
                                                  sep='\t', index=False)
