#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility helpers for the CS-and-rankings data-exploration notebooks
(IPF, ALS-NYGC, ...).

`find_drug_in_run` answers one question: is a given drug present in a loaded
connectivity-score run, and where does it rank?
"""
import re

import numpy as np
import pandas as pd

_TIMEPOINTS = ('6h_24h', '24h', '6h')


def _normalize(name):
    """Lower-case + strip every non-alphanumeric char (the same rule `load_cs`
    uses to build `drug_key`)."""
    return re.sub(r'[^a-z0-9]', '', str(name).lower())


def _split_pert_id(pid):
    """Split '<drug>_<timepoint>' into (drug, perturbation_time)."""
    for tp in _TIMEPOINTS:
        if str(pid).endswith('_' + tp):
            return pid[:-(len(tp) + 1)], tp
    return pid, np.nan


def find_drug_in_run(data, drug_name, perturbation_time=None, contains=False,
                     verbose=True):
    """Inspect whether `drug_name` is present in a CS run, and where it ranks.

    Matching is done on a normalized drug key (lower-cased, non-alphanumeric
    stripped) -- the same normalization `load_cs` uses for `drug_key` -- so the
    query's casing and punctuation are irrelevant: 'Edaravone', 'edaravone' and
    'ICI-199441' all match. If the run frame doesn't already carry the derived
    `drug` / `drug_key` columns, they are reconstructed from `pert_id`, so this
    also works on a raw CS .tsv loaded via loader.load_cs_run().

    Parameters
    ----------
    data : pd.DataFrame
        A CS run. Needs a `connectivity_score` column and either the derived
        `drug` / `drug_key` columns (as added by the notebooks' `load_cs`) or a
        `pert_id` column to rebuild them from.
    drug_name : str
        Drug name to look for; normalized the same way as `drug_key`.
    perturbation_time : str | None
        Optional timepoint filter ('6h', '24h', '6h_24h').
    contains : bool
        False (default) -> exact key match; True -> substring match (handy for
        partial names / salt forms).
    verbose : bool
        Print a one-line found / not-found summary.

    Returns
    -------
    pd.DataFrame
        Matching rows, most-negative connectivity_score first, with a `rank`
        column = the signature's 1-based rank within the whole run by
        connectivity_score (rank 1 = most-negative = strongest reverser, per the
        NetCoS sign convention). Empty frame if the drug is absent.
    """
    df = data.copy()
    if 'connectivity_score' not in df.columns:
        raise ValueError("run frame has no 'connectivity_score' column")

    # Rebuild the derived columns if this is a raw CS frame.
    if 'drug' not in df.columns or 'perturbation_time' not in df.columns:
        if 'pert_id' not in df.columns:
            raise ValueError("run frame needs either 'drug'/'perturbation_time' "
                             "or a 'pert_id' column to look drugs up")
        split = df['pert_id'].apply(_split_pert_id)
        df['drug'] = split.apply(lambda x: x[0])
        df['perturbation_time'] = split.apply(lambda x: x[1])
    if 'drug_key' not in df.columns:
        df['drug_key'] = df['drug'].map(_normalize)

    df['rank'] = df['connectivity_score'].rank(method='min').astype(int)

    key = _normalize(drug_name)
    if contains:
        hit = df[df['drug_key'].str.contains(key, na=False)]
    else:
        hit = df[df['drug_key'] == key]
    if perturbation_time is not None:
        hit = hit[hit['perturbation_time'] == perturbation_time]

    cols = [c for c in ['drug', 'perturbation_time', 'connectivity_score',
                        'cs_p_value', 'rank', 'pearson', 'spearman', 'cos_sim']
            if c in hit.columns]
    hit = hit[cols].sort_values('connectivity_score')

    if verbose:
        where = '' if perturbation_time is None else f' at {perturbation_time}'
        if hit.empty:
            print(f"'{drug_name}' (key='{key}') not in run{where}.")
        else:
            print(f"'{drug_name}' found: {len(hit)} signature(s){where} "
                  f"(run has {len(df)} rows).")
    return hit
