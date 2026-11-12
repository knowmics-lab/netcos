#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUROC of approved vs negative asthma drugs (repoDB pairs with a curated
trial stop reason) for the per-drug scores of 03_srges_per_drug.py and for
the drug-collapsed (TSR/LMM) connectivity-score runs.

Design
------
- Labels: repodb_signature_lincs_pairs_stop_reason.tsv rows with
  ind_id == ASTHMA_CUI. Positives = 'approved'. Negatives = 'not_approved'
  pairs whose pair_stop_class is in NEG_CLASSES (a trial stopped for
  efficacy, safety, or both). Two negative sets: 'stop_efficacy_safety'
  (pair_stop_class) and 'stop_efficacy_safety_enrolled'
  (clinical_stop_enrolled, i.e. Terminated or Suspended trials only, so
  that participants were dosed). Every other not-approved pair (non-clinical,
  interim_unspecified or unknown stop) and every 'conflict' pair is excluded.
- Matching: per-instance route through the L1000 compound ids
  (lincs_pert_ids, '|'-separated) to their L1000 name, the compound key of
  03_srges_per_drug.py; drug-collapsed route on the normalised L1000 name
  (lincs_pert_inames, lower case, alphanumerics only; alias
  beclometasone -> beclomethasonedipropionate), pooled 6h_24h rows.
  A drug matched to several compounds takes the mean score.
- Drug sets: 'all_scored' = every labelled drug with a score in that
  route/configuration; 'matched_to_tsr' (per-instance route only) = the
  labelled drugs that also have a score in the TSR run of the same
  configuration, so that both routes are compared on the same drugs.
- Score = -(per-drug score): AUROC > 0.5 means approved drugs have more
  negative (more reversing) scores than the negatives. AUROC = Mann-Whitney U / (n_pos n_neg),
  exact two-sided p-value; 95% percentile bootstrap (2,000 resamples within
  each class, seed 0).
- Drug-collapsed runs: the most recent cs_batch run of conf.DISEASE on drug
  input 'all' per configuration, via loader.load_cs_run.

Inputs
------
- P.oout_dir/<slug>_srges_per_drug.tsv
- LABELS_FILE
- conf.cs_log_filename (drug-collapsed runs)

Outputs
-------
- P.oout_dir/<slug>_repodb_auroc_srges.tsv
    one row per (negative_set, drug_set, route, config, collapse): n_approved,
    n_negative, AUROC, CI95_low, CI95_high, p_MannWhitney.
- P.oout_dir/<slug>_repodb_labelled_drug_scores.tsv
    every approved or negative drug (label, pair_stop_class,
    clinical_stop_enrolled) with its score in every (route, config, collapse)
    and its number of instances / cell lines in the per-instance route.

Usage
-----
    cd src && PYTHONPATH=. python case_studies/asthma_GSE152004/04_repodb_auroc.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

from case_conf import P
import conf
from loader import load_cs_run

ASTHMA_CUI = 'C0004096'
LABELS_FILE = conf.OOUT_DIR / 'validations' / 'drug_disease_pairs_validation' / 'repodb_signature_lincs_pairs_stop_reason.tsv'
NEG_CLASSES = ['efficacy', 'safety', 'efficacy_and_safety']
NEGATIVE_SETS = {'stop_efficacy_safety': 'pair_stop_class', 'stop_efficacy_safety_enrolled': 'clinical_stop_enrolled'}
NAME_ALIAS = {'beclometasone': 'beclomethasonedipropionate'}
TSR_CONFIGS = {'DEG': (0, 0, 0), 'DEG_LM': (0, 1, 0), 'mith': (1, 0, 0), 'mith_LM': (1, 1, 0), 'mith_pw': (1, 0, 1)}


def key(s):
    return pd.Series(s).astype(str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True).values


def auroc(d):
    pos = -d.loc[d['label'] == 'approved', 'score'].values
    neg = -d.loc[d['label'] == 'negative', 'score'].values
    if len(pos) == 0 or len(neg) == 0:
        return dict(n_approved=len(pos), n_negative=len(neg))
    mw = stats.mannwhitneyu(pos, neg, alternative='two-sided', method='exact')
    rng = np.random.default_rng(0)
    boot = [stats.mannwhitneyu(rng.choice(pos, len(pos)), rng.choice(neg, len(neg))).statistic
            / (len(pos) * len(neg)) for _ in range(2000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(n_approved=len(pos), n_negative=len(neg), AUROC=mw.statistic / (len(pos) * len(neg)),
                CI95_low=lo, CI95_high=hi, p_MannWhitney=mw.pvalue)


def labels():
    lab = pd.read_csv(LABELS_FILE, sep='\t')
    lab = lab[lab['ind_id'] == ASTHMA_CUI].copy()
    neg = (lab['label'] == 'not_approved') & lab['pair_stop_class'].isin(NEG_CLASSES)
    print(f'asthma pairs by repoDB label: {lab["label"].value_counts().to_dict()}')
    print('not-approved pairs by pair_stop_class:',
          lab.loc[lab['label'] == 'not_approved', 'pair_stop_class'].value_counts(dropna=False).to_dict())
    lab = lab[(lab['label'] == 'approved') | neg].copy()
    lab.loc[lab['label'] == 'not_approved', 'label'] = 'negative'
    cols = ['drugbank_id', 'drug_name', 'label', 'pair_stop_class', 'clinical_stop_enrolled',
            'lincs_pert_ids', 'lincs_pert_inames']
    return lab[cols].reset_index(drop=True)


def per_instance_scores(lab):
    s = pd.read_csv(P.oout_dir / f'{P.slug}_srges_per_drug.tsv', sep='\t')
    names = pd.read_csv(conf.LINCS_METADATA_PATH, usecols=['pert_id', 'pert_iname']).drop_duplicates('pert_id')
    ex = lab.assign(pert_id=lab['lincs_pert_ids'].astype(str).str.split('|')).explode('pert_id')
    ex = ex.merge(names, on='pert_id', how='left').drop_duplicates(['drugbank_id', 'drug_name', 'label', 'pert_iname'])
    m = ex.merge(s, on='pert_iname', how='left')
    g = m.groupby(['drugbank_id', 'drug_name', 'label', 'pair_stop_class', 'clinical_stop_enrolled',
                   'config', 'collapse'], dropna=False)
    out = g.agg(score=('score', 'mean'), n_instances=('n_instances', 'sum'),
                n_cell_lines=('n_cell_lines', 'max'), n_names=('pert_iname', 'nunique')).reset_index()
    out = out[out['config'].notna()]
    out.insert(0, 'route', 'per_instance_all_lines')
    return out


def tsr_scores(lab):
    ex = lab.assign(name=lab['lincs_pert_inames'].astype(str).str.split('|')).explode('name')
    ex['k'] = key(ex['name'])
    ex['k'] = ex['k'].replace(NAME_ALIAS)
    rows = []
    for label, (m, lm, pw) in TSR_CONFIGS.items():
        try:
            df, path, run_id = load_cs_run(disease=conf.DISEASE, landmark_disease=False, drug_run_id='all', mith=m,
                                           cs_on_LM=lm, CS_ON_PATHWAYS=pw, CS_METHOD=conf.CS_METHOD, verbose=False)
        except FileNotFoundError:
            continue
        pooled = df[df['pert_id'].astype(str).str.endswith('_6h_24h')].copy()
        pooled['k'] = key(pooled['pert_id'].str[:-len('_6h_24h')])
        sc = pooled.groupby('k')['connectivity_score'].mean()
        r = ex.assign(score=ex['k'].map(sc)).groupby(['drugbank_id', 'drug_name', 'label', 'pair_stop_class',
                                                      'clinical_stop_enrolled'], dropna=False)['score'].mean().reset_index()
        r['config'], r['collapse'], r['cs_run_id'] = label, 'tsr_lmm_6h_24h', run_id
        rows.append(r)
    out = pd.concat(rows, ignore_index=True)
    out.insert(0, 'route', 'tsr_lmm_collapsed')
    return out


if __name__ == '__main__':
    lab = labels()
    print(f'labelled asthma drugs kept: {lab["label"].value_counts().to_dict()}; negatives: '
          + ', '.join(f'{r.drug_name} ({r.pair_stop_class}, enrolled: {r.clinical_stop_enrolled})'
                      for r in lab[lab.label == 'negative'].itertuples()))
    scores = pd.concat([per_instance_scores(lab), tsr_scores(lab)], ignore_index=True)
    res = []
    tsr_ids = {(c, k) for c, k in zip(scores.loc[scores.route == 'tsr_lmm_collapsed'].dropna(subset=['score'])['config'],
                                         scores.loc[scores.route == 'tsr_lmm_collapsed'].dropna(subset=['score'])['drugbank_id'])}
    for nset, col in NEGATIVE_SETS.items():
        keep = (scores['label'] == 'approved') | scores[col].isin(NEG_CLASSES)
        for (route, config, coll), d in scores[keep].groupby(['route', 'config', 'collapse']):
            d = d.dropna(subset=['score'])
            base = {'negative_set': nset, 'route': route, 'config': config, 'collapse': coll}
            res.append({**base, 'drug_set': 'all_scored', **auroc(d)})
            if route == 'per_instance_all_lines':
                dm = d[[(config, i) in tsr_ids for i in d['drugbank_id']]]
                res.append({**base, 'drug_set': 'matched_to_tsr', **auroc(dm)})
    res = pd.DataFrame(res)
    res.to_csv(P.oout_dir / f'{P.slug}_repodb_auroc_srges.tsv', sep='\t', index=False)
    scores.to_csv(P.oout_dir / f'{P.slug}_repodb_labelled_drug_scores.tsv', sep='\t', index=False)
    print(res.round(3).to_string(index=False))
