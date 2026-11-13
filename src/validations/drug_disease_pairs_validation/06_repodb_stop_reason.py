#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stop reason of the failed trials behind the not-approved repoDB pairs, and the
L-DDDEP labels (LINCS Drug Disease Differential Expression Pairs, working name)
of every drug-disease pair of the 04_ table.

Design
------
- 'DetailedStatus' of repoDB v2 (conf.REPODB_FILE) is the ClinicalTrials.gov
  "Why Study Stopped" text of each failed (Terminated / Withdrawn / Suspended)
  record. Texts are compared after collapsing whitespace.
- Stop reason of a text, two curated inputs in this slug's output folder:
  1. WHY_STOPPED_FILE (why_stopped, stop_reason, note): one category per
     distinct text. The texts of the not-approved pairs at >= 10 Psamples come
     from one LLM-assisted pass cross-checked against keyword rules
     (disagreements resolved by hand); the texts added for the pairs at 5-9
     Psamples are a single pass, not reviewed (see 'note'). Categories:
     efficacy, safety, efficacy_and_safety, interim_unspecified, recruitment,
     business_funding, regulatory_admin, covid, other_evidence,
     indirect_evidence, completed_not_failure, unclear.
  2. REVIEW_FILE (why_stopped, ..., lorenzo): Lorenzo's decisions on the
     doubtful texts (other_evidence and stops that read as positive). The
     first word of 'lorenzo' is the decision and overrides the category:
     uninformative, completed_not_failure, negative_efficacy_external
     (-> indirect_evidence, not a negative), negative_safety_external
     (-> safety_external, a negative).
  Records without text are 'missing'; texts in neither file 'not_curated'.
- Trial class (trial_class) from the final stop reason:
  efficacy / safety / efficacy_and_safety / safety_external -> negative_<reason>;
  indirect_evidence and completed_not_failure kept as such; every other
  category -> uninformative; missing and not_curated kept as such.
- Psamples of a pair (psamples): the smaller arm (case or control) of the
  largest signature dataset of the disease, max over DiSignAtlas and CREEDS
  (columns dsa_/creeds_biggest_ctrl_case of the 04_ table; 'ctrl|case', one
  pair of numbers per tied dataset). max_psamples of 04_ is its threshold bin.
- L-DDDEP label of a pair (ldddep_label): 'positive' for repoDB approved,
  'conflict' for conflict pairs, and for not-approved pairs 'negative' if at
  least one of its trials (any status) has a negative_* trial class, else
  'excluded'. ldddep_negative_reasons lists the negative trial classes.
- Pair stop classes of the previous version, kept unchanged (computed on the
  category of WHY_STOPPED_FILE before the review; trial_stop_reasons too):
  pair_stop_class = the first of efficacy_and_safety, efficacy, safety,
  interim_unspecified among the pair's trials, else 'non_clinical' if at
  least one trial has a non-clinical reason, else 'unknown'; empty for
  approved pairs. clinical_stop_enrolled = the same on Terminated and
  Suspended trials only. all_withdrawn = every failed trial is Withdrawn.

Inputs
------
- conf.REPODB_FILE.
- OUT_DIR/repodb_signature_lincs_pairs.tsv (from 04_).
- OUT_DIR/repodb_why_stopped_classification.tsv (curated by hand).
- OUT_DIR/review_other_evidence_and_positive_stops.tsv (Lorenzo's decisions).

Outputs (OUT_DIR = OOUT_DIR/validations/drug_disease_pairs_validation/)
----------------------------------------------------------------------
- ldddep_pairs_by_trial.tsv : THE L-DDDEP DATASET. One row per (pair, failed
  trial); approved pairs have one row with empty trial columns. Pair columns
  (drug, disease, psamples, repoDB label, ldddep_label, LINCS mapping) and
  trial columns (NCT, status, phase, why_stopped full text, stop_reason,
  label_source, lorenzo_decision, trial_class).
- ldddep_summary.tsv : labelled pairs, diseases and drugs at >= 5/10/20/50 Psamples.
- repodb_stop_reason_by_trial.tsv, repodb_signature_lincs_pairs_stop_reason.tsv,
  repodb_stop_reason_summary_ge10Psamples.tsv : as in the previous version
  (the pair table also gains psamples, ldddep_label, ldddep_negative_reasons).

Usage
-----
    PYTHONPATH=. python validations/drug_disease_pairs_validation/06_repodb_stop_reason.py   (from src/)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SRC_DIR))

from conf import REPODB_FILE, OOUT_DIR

SLUG = "drug_disease_pairs_validation"
OUT_DIR = OOUT_DIR / "validations" / SLUG
WHY_STOPPED_FILE = OUT_DIR / "repodb_why_stopped_classification.tsv"
REVIEW_FILE = OUT_DIR / "review_other_evidence_and_positive_stops.tsv"
PSAMPLES = 10
THRESHOLDS = [5, 10, 20, 50]
CLINICAL = ['efficacy_and_safety', 'efficacy', 'safety', 'interim_unspecified']
NON_CLINICAL = {'recruitment', 'business_funding', 'regulatory_admin', 'covid', 'other_evidence'}
REVIEW_DECISION = {'uninformative': 'uninformative', 'completed_not_failure': 'completed_not_failure',
                   'negative_efficacy_external': 'indirect_evidence', 'negative_safety_external': 'safety_external'}
NEGATIVE = {'efficacy': 'negative_efficacy', 'safety': 'negative_safety',
            'efficacy_and_safety': 'negative_efficacy_and_safety', 'safety_external': 'negative_safety_external'}
KEEP = {'indirect_evidence', 'completed_not_failure', 'missing', 'not_curated'}
PAIR_COLS = ['drugbank_id', 'drug_name', 'ind_id', 'ind_name', 'psamples', 'max_psamples', 'label',
             'ldddep_label', 'ldddep_negative_reasons', 'n_failed_trials', 'signature_databases',
             'lincs_pert_inames', 'lincs_n_signatures', 'lincs_n_cell_lines']
TRIAL_COLS = ['NCT', 'status', 'phase', 'why_stopped', 'stop_reason', 'label_source', 'lorenzo_decision', 'trial_class']


def norm(s):
    return s.astype('string').str.replace(r'\s+', ' ', regex=True).str.strip()


def pair_class(reasons):
    for c in CLINICAL:
        if c in reasons:
            return c
    return 'non_clinical' if NON_CLINICAL & set(reasons) else 'unknown'


def smaller_arm(cell):
    if pd.isna(cell):
        return np.nan
    v = [int(x) for x in str(cell).split('|')]
    assert len(v) % 2 == 0, f'odd number of sample sizes: {cell}'
    return max(min(v[i], v[i + 1]) for i in range(0, len(v) - 1, 2))


def trial_class(reason):
    if reason in NEGATIVE:
        return NEGATIVE[reason]
    return reason if reason in KEEP else 'uninformative'


if __name__ == '__main__':
    # ---- stop reason per text: curated table, then Lorenzo's review on top
    dec = pd.read_csv(WHY_STOPPED_FILE, sep='\t', keep_default_na=False)
    dec['why_stopped'] = norm(dec['why_stopped'])
    assert not dec['why_stopped'].duplicated().any(), 'duplicated texts in the decision table'
    rev = pd.read_csv(REVIEW_FILE, sep='\t', keep_default_na=False)
    rev['why_stopped'] = norm(rev['why_stopped'])
    rev['decision'] = rev['lorenzo'].str.split(r'[\s:(]', regex=True).str[0]
    bad = set(rev['decision']) - set(REVIEW_DECISION)
    assert not bad, f'unknown review decisions: {bad}'
    reason = dict(zip(dec['why_stopped'], dec['stop_reason']))
    review = dict(zip(rev['why_stopped'], rev['decision'].map(REVIEW_DECISION)))
    review_text = dict(zip(rev['why_stopped'], rev['lorenzo']))

    repo = pd.read_csv(REPODB_FILE)
    failed = repo[repo['status'] != 'Approved'].copy()
    failed['why_stopped'] = norm(failed['DetailedStatus'])
    src, rsn = [], []
    for t in failed['why_stopped']:
        if pd.isna(t) or t == '':
            src.append('missing'); rsn.append('missing')
        elif t in review:
            src.append('lorenzo_review'); rsn.append(review[t])
        elif t in reason:
            src.append('curated_table'); rsn.append(reason[t])
        else:
            src.append('not_curated'); rsn.append('not_curated')
    failed['label_source'] = src
    failed['stop_reason'] = rsn
    # category of the curated table before the review, for pair_stop_class (as in the previous version)
    failed['stop_reason_table'] = [reason.get(t, 'missing' if pd.isna(t) or t == '' else 'not_curated')
                                   for t in failed['why_stopped']]
    failed['lorenzo_decision'] = failed['why_stopped'].map(review_text).fillna('')
    failed['trial_class'] = failed['stop_reason'].map(trial_class)
    # one row per NCT: every record of an NCT must carry the same status and text
    per_nct = failed.groupby('NCT').agg(n_status=('status', 'nunique'), n_text=('why_stopped', 'nunique'))
    assert (per_nct[['n_status', 'n_text']] <= 1).all().all(), 'an NCT has records with different status or text'
    trial_info = (failed.drop_duplicates('NCT')
                  .set_index('NCT')[['status', 'phase', 'why_stopped', 'stop_reason', 'label_source',
                                     'lorenzo_decision', 'trial_class', 'stop_reason_table']])

    # ---- pairs of 04_
    pairs = pd.read_csv(OUT_DIR / 'repodb_signature_lincs_pairs.tsv', sep='\t')
    pairs['psamples'] = pairs[['dsa_biggest_ctrl_case', 'creeds_biggest_ctrl_case']].map(smaller_arm).max(axis=1)
    rows, per_pair = [], []
    for _, r in pairs.iterrows():
        ncts = [] if pd.isna(r['failed_trials']) else str(r['failed_trials']).split('|')
        info = trial_info.reindex(ncts)
        info['stop_reason'] = info['stop_reason'].fillna('not_curated')
        info['trial_class'] = info['trial_class'].fillna('not_curated')
        info['stop_reason_table'] = info['stop_reason_table'].fillna('not_curated')
        for nct, t in info.iterrows():
            rows.append({'drugbank_id': r['drugbank_id'], 'ind_id': r['ind_id'], 'NCT': nct,
                         **t.drop('stop_reason_table').to_dict()})
        if r['label'] == 'approved' or not ncts:
            lab = {'approved': 'positive', 'conflict': 'conflict'}.get(r['label'], 'excluded')
            per_pair.append({'n_failed_trials': len(ncts), 'trial_stop_reasons': '', 'pair_stop_class': '',
                             'clinical_stop_enrolled': '', 'all_withdrawn': '',
                             'ldddep_label': lab, 'ldddep_negative_reasons': ''})
            continue
        reasons = info['stop_reason_table'].tolist()
        enrolled = info.loc[info['status'].isin(['Terminated', 'Suspended']), 'stop_reason_table'].tolist()
        negs = sorted({c for c in info['trial_class'] if c.startswith('negative_')})
        lab = 'conflict' if r['label'] == 'conflict' else ('negative' if negs else 'excluded')
        per_pair.append({'n_failed_trials': len(ncts), 'trial_stop_reasons': '|'.join(reasons),
                         'pair_stop_class': pair_class(reasons),
                         'clinical_stop_enrolled': pair_class(enrolled) if enrolled else 'none_enrolled',
                         'all_withdrawn': bool((info['status'] == 'Withdrawn').all()),
                         'ldddep_label': lab, 'ldddep_negative_reasons': '|'.join(negs)})
    out = pd.concat([pairs, pd.DataFrame(per_pair, index=pairs.index)], axis=1)
    by_trial = pd.DataFrame(rows, columns=['drugbank_id', 'ind_id'] + TRIAL_COLS)

    # ---- the L-DDDEP dataset: pair x trial, approved pairs with empty trial columns
    ldddep = out[PAIR_COLS].merge(by_trial, on=['drugbank_id', 'ind_id'], how='left')
    ldddep = ldddep.sort_values(['ind_name', 'ldddep_label', 'drug_name', 'NCT'], na_position='last')

    # ---- summary per Psamples threshold
    summ = []
    for th in THRESHOLDS:
        q = out[out['psamples'] >= th]
        lab = q[q['ldddep_label'].isin(['positive', 'negative'])]
        per_dis = lab.pivot_table(index='ind_id', columns='ldddep_label', values='drugbank_id',
                                  aggfunc='nunique', fill_value=0).reindex(columns=['positive', 'negative'], fill_value=0)
        summ.append({'min_psamples': th, **q['ldddep_label'].value_counts().to_dict(),
                     'diseases': q['ind_id'].nunique(),
                     'diseases_with_negative': int((per_dis['negative'] > 0).sum()),
                     'diseases_with_both': int(((per_dis['positive'] > 0) & (per_dis['negative'] > 0)).sum()),
                     'diseases_with_ge5_each': int(((per_dis['positive'] >= 5) & (per_dis['negative'] >= 5)).sum()),
                     'drugs_labelled': lab['drugbank_id'].nunique()})
    summ = pd.DataFrame(summ)

    # ---- write
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ldddep.to_csv(OUT_DIR / 'ldddep_pairs_by_trial.tsv', sep='\t', index=False)
    summ.to_csv(OUT_DIR / 'ldddep_summary.tsv', sep='\t', index=False)
    old_trial = out[['drugbank_id', 'drug_name', 'ind_id', 'ind_name', 'label']].merge(by_trial, on=['drugbank_id', 'ind_id'])
    old_trial = old_trial[old_trial['label'] != 'approved']
    old_trial.to_csv(OUT_DIR / 'repodb_stop_reason_by_trial.tsv', sep='\t', index=False)
    out.to_csv(OUT_DIR / 'repodb_signature_lincs_pairs_stop_reason.tsv', sep='\t', index=False)

    na = out[(out['label'] == 'not_approved') & (out['max_psamples'] >= PSAMPLES)]
    summary = pd.concat([
        na['pair_stop_class'].value_counts().rename_axis('class').reset_index(name='n_pairs').assign(view='pair_stop_class'),
        na['clinical_stop_enrolled'].value_counts().rename_axis('class').reset_index(name='n_pairs').assign(view='clinical_stop_enrolled'),
        pd.DataFrame([{'view': 'all_withdrawn', 'class': 'True', 'n_pairs': int((na['all_withdrawn'] == True).sum())}]),
        na['ldddep_label'].value_counts().rename_axis('class').reset_index(name='n_pairs').assign(view='ldddep_label'),
    ])[['view', 'class', 'n_pairs']]
    summary.to_csv(OUT_DIR / f'repodb_stop_reason_summary_ge{PSAMPLES}Psamples.tsv', sep='\t', index=False)
    print(f'not-approved pairs at >= {PSAMPLES} Psamples: {len(na)}')
    print(summary.to_string(index=False))
    print('\nL-DDDEP labels by Psamples threshold:')
    print(summ.to_string(index=False))
    print('\nlabel source of the trials of not-approved pairs at >= 5 Psamples:')
    q5 = ldddep[(ldddep['label'] == 'not_approved') & (ldddep['psamples'] >= 5)].drop_duplicates('NCT')
    print(q5['label_source'].value_counts().to_string())
