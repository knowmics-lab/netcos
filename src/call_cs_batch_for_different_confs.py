# -*- coding: utf-8 -*-
"""
Created on Fri May  8 2026

@author: L-F-S

Loop cs_batch.run_cs_batch_for_conf over a grid of hyperparameters.

- Drops combos that violate the conf.validate_hyperparameters constraints
  (CS_ON_PATHWAYS=True requires cs_mith=1 and cs_on_LM=0).
- Reads logs/cs_runs.tsv. If any combo in the grid is already logged there,
  prompts once whether to skip them; otherwise runs everything without
  further interaction.

Run from src/:
    python call_cs_batch_for_different_confs.py
"""

import itertools
import sys

import pandas as pd

import conf
from conf import cs_log_filename, diseases_of, validate_hyperparameters
from cs_batch import run_cs_batch_for_conf


# ===================== GRID (edit here) =====================
# Single source of truth for the sweep. Mirrors the conf.py knobs.
cell_lines        = ['HEPG2','MCF7','HT29']                                
landmark_diseases = [False, True]
landmark_drugs    = [True]                                    # BinChen drug data is already LM only
CS_METHODs        = ['bin_chen_disease_sorted', 'bin_chen']
cs_miths          = [0, 1]
cs_on_LMs         = [0, 1]
CS_ON_PATHWAYSs   = [False, True]
# ============================================================


def is_compatible(cs_mith, cs_on_LM, CS_ON_PATHWAYS):
    """Wrap conf.validate_hyperparameters in a bool predicate."""
    try:
        validate_hyperparameters(cs_mith, cs_on_LM, CS_ON_PATHWAYS)
        return True
    except ValueError:
        return False


def _to_bool(v):
    """Coerce bool/int/float/str values from cs_runs.tsv into a Python bool."""
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return bool(int(v))
    s = str(v).strip().lower()
    return s in ('true', '1', '1.0')


def already_logged(log_df, cell_line, landmark_disease, landmark_drug,
                   CS_METHOD, cs_mith, cs_on_LM, CS_ON_PATHWAYS):
    """
    True if a row with the same (disease_run_id, drug_run_id, mith, cs_on_LM,
    CS_ON_PATHWAYS, CS_METHOD) tuple already exists in cs_runs.tsv.

    CS_ON_PATHWAYS has historically been written as bool, int, or str depending
    on dtype inference at concat time, so we normalize both sides via _to_bool.
    """
    if log_df.empty:
        return False
    disease_run_id = diseases_of[cell_line] + ('_LM' if landmark_disease else '')
    drug_run_id    = cell_line              + ('_LM' if landmark_drug    else '')
    log_pw_bool = log_df['CS_ON_PATHWAYS'].apply(_to_bool)
    match = (
        (log_df['disease_run_id'] == disease_run_id) &
        (log_df['drug_run_id']    == drug_run_id) &
        (log_df['mith'].astype(int)        == int(cs_mith)) &
        (log_df['cs_on_LM'].astype(int)    == int(cs_on_LM)) &
        (log_pw_bool                       == bool(CS_ON_PATHWAYS)) &
        (log_df['CS_METHOD']               == CS_METHOD)
    )
    return bool(match.any())


def format_combo(combo):
    cell_line, ld_disease, ld_drug, method, mith, on_lm, on_pw = combo
    return (f"cell_line={cell_line} landmark_disease={ld_disease} "
            f"landmark_drug={ld_drug} CS_METHOD={method} "
            f"cs_mith={mith} cs_on_LM={on_lm} CS_ON_PATHWAYS={on_pw}")


def main():
    all_combos = list(itertools.product(
        cell_lines, landmark_diseases, landmark_drugs,
        CS_METHODs, cs_miths, cs_on_LMs, CS_ON_PATHWAYSs))

    valid   = [c for c in all_combos if is_compatible(c[4], c[5], c[6])]
    invalid = [c for c in all_combos if c not in valid]

    log_df = (pd.read_csv(cs_log_filename, sep='\t')
              if cs_log_filename.exists() and cs_log_filename.stat().st_size > 0
              else pd.DataFrame())

    duplicates = [c for c in valid if already_logged(log_df, *c)]
    fresh      = [c for c in valid if c not in duplicates]

    print(f"Grid total combos       : {len(all_combos)}")
    print(f"  invalid (skipped)     : {len(invalid)}")
    print(f"  valid                 : {len(valid)}")
    print(f"    already in cs_runs  : {len(duplicates)}")
    print(f"    fresh               : {len(fresh)}")

    to_run = fresh
    if duplicates:
        print("\nThe following valid combos already have a row in cs_runs.tsv:")
        for c in duplicates:
            print(f"  - {format_combo(c)}")
        ans = input("\nSkip these already-logged combos? [y/N]: ").strip().lower()
        if ans != 'y':
            to_run = fresh + duplicates

    if not to_run:
        print("\nNothing to run. Exiting.")
        return


    for i, combo in enumerate(to_run, 1):
        cell_line, ld_disease, ld_drug, method, mith, on_lm, on_pw = combo
        print(f"\n=== [{i}/{len(to_run)}] {format_combo(combo)} ===")
        try:
            run_cs_batch_for_conf(
                cell_line=cell_line,
                landmark_disease=ld_disease,
                landmark_drug=ld_drug,
                CS_METHOD=method,
                cs_mith=mith,
                cs_on_LM=on_lm,
                CS_ON_PATHWAYS=on_pw,
            )
        except Exception as e:
            print(f"!!! Combo failed ({format_combo(combo)}): {e}", file=sys.stderr)
            ans = input("Continue with remaining combos? [y/N]: ").strip().lower()
            if ans != 'y':
                raise

    print("\nAll combos done.")


if __name__ == '__main__':
    main()
