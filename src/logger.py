# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:10:01 2026

@author: los4
"""
import fcntl
import pandas as pd
from pathlib import Path

import conf


def append_run_metadata(metadata_path, row_dict):
    """
    Append one row to a run log TSV, creating the file if needed.

    The read-concat-write runs under an exclusive fcntl lock on
    conf.LOGS_DIR / 'cs_runs.tsv.lock', so concurrent writers never lose rows.
    """
    metadata_path = Path(metadata_path)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = Path(conf.LOGS_DIR) / 'cs_runs.tsv.lock'
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    row_df = pd.DataFrame([row_dict])

    with open(lock_path, 'a') as lock_fh:
        fcntl.flock(lock_fh, fcntl.LOCK_EX)
        try:
            if metadata_path.exists() and metadata_path.stat().st_size > 0:
                old_df = pd.read_csv(metadata_path, sep="\t")
                out_df = pd.concat([old_df, row_df], ignore_index=True)
            else:
                out_df = row_df
            out_df.to_csv(metadata_path, sep="\t", index=False)
        finally:
            fcntl.flock(lock_fh, fcntl.LOCK_UN)
