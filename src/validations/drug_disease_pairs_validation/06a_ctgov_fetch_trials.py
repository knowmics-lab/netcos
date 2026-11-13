#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download the ClinicalTrials.gov record of every failed trial behind the repoDB
pairs of the 04_ table, as raw JSON, for the L-DDDEP labels built by 06_.

Design
------
- Why: repoDB pairs every intervention of a failed trial with the trial's
  condition, whatever the arm type. Comparators, positive controls and
  background therapies become not-approved pairs (e.g. NCT01721135: a thorough
  QT study of GSK2190915 for asthma, with moxifloxacin as positive control,
  becomes "moxifloxacin - asthma, not approved"). The arm type of each
  intervention, the enrolment and the stop reason are in the trial record.
- Which trials: every NCT listed in 'failed_trials' of
  repodb_signature_lincs_pairs.tsv (output of 04_), any label.
- Source: ClinicalTrials.gov API v2, https://clinicaltrials.gov/api/v2/studies/<NCT>
  (one JSON per study). Saved unchanged as DATA_DIR/ClinicalTrials/studies/<NCT>.json
  (raw downloaded data). Files already on disk are not downloaded again, so the
  script can be re-run after an interruption.
- Requests: one at a time, PAUSE_SEC between them, RETRIES attempts with
  back-off on errors and on HTTP 429 / 5xx; HTTP 404 is recorded and not retried.
- Checks at the start: number of pairs and NCTs, NCT format, NCTs already on
  disk. At the end: a log of the run and a check that every saved file parses
  and holds the requested NCT.

Inputs
------
- OOUT_DIR/validations/drug_disease_pairs_validation/repodb_signature_lincs_pairs.tsv (04_).

Outputs
-------
- DATA_DIR/ClinicalTrials/studies/<NCT>.json : raw API records.
- OOUT_DIR/validations/drug_disease_pairs_validation/ctgov_fetch_log.tsv :
  one row per NCT (nct, status: on_disk / downloaded / not_found / failed,
  http_code, message), rewritten at every run.

Usage
-----
    PYTHONPATH=. python validations/drug_disease_pairs_validation/06a_ctgov_fetch_trials.py   (from src/)
Needs network access to clinicaltrials.gov (a VPN may get HTTP 403).
"""

import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

SRC_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SRC_DIR))

from conf import DATA_DIR, OOUT_DIR

SLUG = "drug_disease_pairs_validation"
OUT_DIR = OOUT_DIR / "validations" / SLUG
CTGOV_DIR = DATA_DIR / "ClinicalTrials" / "studies"
API = "https://clinicaltrials.gov/api/v2/studies/{}"
PAUSE_SEC = 0.3
RETRIES = 4
TIMEOUT_SEC = 30


def fetch(nct, session):
    """Return (status, http_code, message)."""
    for attempt in range(RETRIES):
        try:
            r = session.get(API.format(nct), timeout=TIMEOUT_SEC)
        except requests.RequestException as e:
            code, msg = None, str(e)[:200]
        else:
            code = r.status_code
            if code == 200:
                rec = r.json()
                got = rec.get('protocolSection', {}).get('identificationModule', {}).get('nctId')
                if got != nct:
                    return 'failed', code, f'record holds {got}'
                (CTGOV_DIR / f'{nct}.json').write_text(json.dumps(rec), encoding='utf-8')
                return 'downloaded', code, ''
            if code == 404:
                return 'not_found', code, ''
            msg = r.text[:200].replace('\n', ' ')
            if code not in (429, 500, 502, 503, 504):
                return 'failed', code, msg
        time.sleep(2 ** attempt * 2)
    return 'failed', code, msg


if __name__ == '__main__':
    # ---- checks on the input
    pairs = pd.read_csv(OUT_DIR / 'repodb_signature_lincs_pairs.tsv', sep='\t')
    ncts = sorted({x for v in pairs['failed_trials'].dropna() for x in str(v).split('|')})
    bad = [x for x in ncts if not re.fullmatch(r'NCT\d{8}', x)]
    print(f'pairs: {len(pairs)} ({pairs["label"].value_counts().to_dict()})')
    print(f'distinct failed trials (NCT): {len(ncts)}; malformed: {bad}')
    assert not bad, 'malformed NCT identifiers'
    CTGOV_DIR.mkdir(parents=True, exist_ok=True)
    todo = [x for x in ncts if not (CTGOV_DIR / f'{x}.json').exists()]
    print(f'already on disk: {len(ncts) - len(todo)}; to download: {len(todo)} -> {CTGOV_DIR}')

    # ---- download
    log = [{'nct': x, 'status': 'on_disk', 'http_code': '', 'message': ''} for x in ncts if x not in todo]
    with requests.Session() as session:
        session.headers['Accept'] = 'application/json'
        for i, nct in enumerate(todo, 1):
            status, code, msg = fetch(nct, session)
            log.append({'nct': nct, 'status': status, 'http_code': code, 'message': msg})
            if i % 50 == 0 or status != 'downloaded':
                print(f'{i}/{len(todo)} {nct} {status} {code or ""} {msg}', flush=True)
            if status == 'failed' and code == 403 and i <= 3:
                sys.exit('HTTP 403 from clinicaltrials.gov: access refused (VPN?). Nothing more requested.')
            time.sleep(PAUSE_SEC)

    log = pd.DataFrame(log).sort_values('nct')
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log.to_csv(OUT_DIR / 'ctgov_fetch_log.tsv', sep='\t', index=False)

    # ---- checks on the output
    unreadable = []
    for nct in ncts:
        f = CTGOV_DIR / f'{nct}.json'
        if f.exists():
            try:
                got = json.loads(f.read_text(encoding='utf-8'))['protocolSection']['identificationModule']['nctId']
                if got != nct:
                    unreadable.append(nct)
            except (ValueError, KeyError):
                unreadable.append(nct)
    print('\nrun summary:', log['status'].value_counts().to_dict())
    print(f'on disk and valid: {sum((CTGOV_DIR / f"{x}.json").exists() for x in ncts) - len(unreadable)} / {len(ncts)}; '
          f'unreadable: {unreadable}')
    print(f'log: {OUT_DIR / "ctgov_fetch_log.tsv"}')
