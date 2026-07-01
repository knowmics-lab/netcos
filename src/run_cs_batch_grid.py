#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run the connectivity score over a grid of pre-MITHrIL / pre-CS hyperparameters
for a single disease.

This is a stripped-down sibling of
validations/chembl/call_cs_batch_for_different_confs.py: it sweeps ONLY the
upstream CS knobs and calls cs_batch.run_cs_batch_for_conf for each combo.
There is no ChEMBL / IC50 / sRGES validation here.

It is meant for the collapsed-drug (LMM / DEG) pipeline (IPF, als_NYGC, ...),
where:
  - the drug side is one already-collapsed signature per drug, living in
    CS_IN_DRUG = connectivity_score/input/drug_signature_2025/<cell_line>
    with cell_line='all'  ->  drug_collapsed_before_cs = True;
  - the disease side is connectivity_score/input/disease_signature_2025/<DISEASE>.

Each combo produces one timestamped CS file in
connectivity_score/output/<DISEASE>[ _LM ]/ and one row in logs/cs_runs.tsv
(both written by run_cs_batch_for_conf).

Edit the GRID block below to choose which hyperparameters to sweep.

Run from src/:
    python run_cs_batch_grid.py                      # disease = als_NYGC (default)
    python run_cs_batch_grid.py --disease ipf
    python run_cs_batch_grid.py --dry-run            # list combos, run nothing
    python run_cs_batch_grid.py --run-all            # don't skip combos already logged
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd

# allow `python run_cs_batch_grid.py` from anywhere
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import conf
from conf import cs_log_filename, validate_hyperparameters
from cs_batch import run_cs_batch_for_conf

# ============================== GRID =================================
# Pre-MITHrIL / pre-CS hyperparameters to sweep. Edit freely.
# Combos that violate conf.validate_hyperparameters are dropped automatically
# (cs_mith=0 forbids pathways; cs_on_LM=1 forbids pathways).
CS_METHODs        = ['bin_chen', 'bin_chen_disease_sorted']
cs_miths          = [0, 1]          # 0 = DEG signatures, 1 = MITHrIL-propagated
cs_on_LMs         = [0, 1]          # 1 = restrict to landmark genes
CS_ON_PATHWAYSs   = [True, False]         # True needs cs_mith=1 & cs_on_LM=0
landmark_diseases = [False]         # adds '_LM' to the disease run name
rank_ons          = ['magnitude']   # {'magnitude', 'p_value'}

# Fixed for the collapsed-drug (LMM/DEG) pipeline. Override --disease on the CLI.
DEFAULT_DISEASE = 'ipf'
CELL_LINE       = 'all'    # drug input dir: drug_signature_2025/all
LANDMARK_DRUG   = False    # the 'all' DEG drug data is not landmark-only
DRUG_COLLAPSED  = True     # signature ids ARE pert_ids (see conf.drug_collapsed_before_cs)
# =====================================================================

# fields, in sweep order, used for labels and the itertools.product
GRID_FIELDS = ['CS_METHOD', 'cs_mith', 'cs_on_LM', 'CS_ON_PATHWAYS',
               'landmark_disease', 'rank_on']


def build_combos():
    """Cartesian product of the grid, dropping conf-invalid combos."""
    all_combos = [dict(zip(GRID_FIELDS, vals)) for vals in itertools.product(
        CS_METHODs, cs_miths, cs_on_LMs, CS_ON_PATHWAYSs,
        landmark_diseases, rank_ons)]
    valid, invalid = [], []
    for c in all_combos:
        try:
            validate_hyperparameters(c['cs_mith'], c['cs_on_LM'], c['CS_ON_PATHWAYS'])
            valid.append(c)
        except ValueError:
            invalid.append(c)
    return valid, invalid


def fmt(disease, c):
    return (f"disease={disease} CS_METHOD={c['CS_METHOD']} cs_mith={c['cs_mith']} "
            f"cs_on_LM={c['cs_on_LM']} CS_ON_PATHWAYS={c['CS_ON_PATHWAYS']} "
            f"landmark_disease={c['landmark_disease']} rank_on={c['rank_on']}")


def _to_bool(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return bool(int(v))
    return str(v).strip().lower() in ('true', '1', '1.0')


def already_logged(log_df, disease, c):
    """True if cs_runs.tsv already has a row for this disease+combo."""
    if log_df.empty:
        return False
    disease_run_id = disease + ('_LM' if c['landmark_disease'] else '')
    drug_run_id    = CELL_LINE + ('_LM' if LANDMARK_DRUG else '')
    match = (
        (log_df['disease_run_id'] == disease_run_id) &
        (log_df['drug_run_id']    == drug_run_id) &
        (log_df['mith'].astype(int)     == int(c['cs_mith'])) &
        (log_df['cs_on_LM'].astype(int) == int(c['cs_on_LM'])) &
        (log_df['CS_ON_PATHWAYS'].apply(_to_bool) == bool(c['CS_ON_PATHWAYS'])) &
        (log_df['CS_METHOD']      == c['CS_METHOD'])
    )
    if 'rank_on' in log_df.columns:
        match &= (log_df['rank_on'] == c['rank_on'])
    return bool(match.any())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--disease', default=DEFAULT_DISEASE,
                    help=f"disease run name (default: {DEFAULT_DISEASE})")
    ap.add_argument('--dry-run', action='store_true',
                    help="list the combos that would run, then exit")
    ap.add_argument('--run-all', action='store_true',
                    help="run every valid combo, even if already in cs_runs.tsv")
    args = ap.parse_args(argv)

    valid, invalid = build_combos()
    
    


    log_df = (pd.read_csv(cs_log_filename, sep='\t')
              if cs_log_filename.exists() and cs_log_filename.stat().st_size > 0
              else pd.DataFrame())

    duplicates = [] if args.run_all else [c for c in valid
                                          if already_logged(log_df, args.disease, c)]
    
    to_run = [c for c in valid if c not in duplicates]

    print("\nINVALID COMBOS:")
    for c in invalid:
        print(f"    - {fmt(args.disease, c)}")
    
    print("\nALREADY LOGGED COMBOS:")
    for c in duplicates:
        print(f"    - {fmt(args.disease, c)}")

    print(f"disease            : {args.disease}")
    print(f"drug input         : drug_signature_2025/{CELL_LINE}"
          f"{'_LM' if LANDMARK_DRUG else ''}  (drug_collapsed={DRUG_COLLAPSED})")
    print(f"grid combos total  : {len(valid) + len(invalid)}")
    print(f"  invalid (skipped): {len(invalid)}")
    print(f"  valid            : {len(valid)}")
    print(f"  already logged   : {len(duplicates)}")
    print(f"  to run           : {len(to_run)}")
    for c in to_run:
        print(f"    - {fmt(args.disease, c)}")

    if args.dry_run or not to_run:
        if not to_run:
            print("\nNothing to run.")
        return

    for i, c in enumerate(to_run, 1):
        print(f"\n=== [{i}/{len(to_run)}] {fmt(args.disease, c)} ===")
        try:
            run_cs_batch_for_conf(
                DISEASE=args.disease,
                cell_line=CELL_LINE,
                landmark_drug=LANDMARK_DRUG,
                drug_collapsed=DRUG_COLLAPSED,
                landmark_disease=c['landmark_disease'],
                CS_METHOD=c['CS_METHOD'],
                cs_mith=c['cs_mith'],
                cs_on_LM=c['cs_on_LM'],
                CS_ON_PATHWAYS=c['CS_ON_PATHWAYS'],
                rank_on=c['rank_on'],
            )
        except Exception as e:
            print(f"!!! combo failed ({fmt(args.disease, c)}): {e}", file=sys.stderr)
            ans = input("Continue with remaining combos? [y/N]: ").strip().lower()
            if ans != 'y':
                raise

    print("\nAll combos done.")


if __name__ == '__main__':
    main()
