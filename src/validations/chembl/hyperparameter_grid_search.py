# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 10:40:00 2026

Cross-validated hyperparameter search for ChEMBL IC50 vs connectivity score validation.

This script is intentionally built on top of final_step_correlations_IC50_CS.py.
It reuses its data-loading, collapsing, classification, and plotting helpers, then:

1. evaluates MANY fixed hyperparameter combinations under out-of-fold test evaluation
2. stores pooled out-of-fold predictions per configuration
3. summarizes precision / recall / F1 / Spearman across configurations
4. plots a matrix of fig3-like scatter panels, one panel per configuration
5. optionally performs nested CV to pick a single best configuration

Important scope note
--------------------
This script only tunes hyperparameters that can be evaluated from already-available
CS / IC50 tables. Hyperparameters that require recomputing the upstream CS itself
(e.g. LM vs all genes, DEG vs MITH, pathway-level vs gene-level signatures, different
CS runs) can still be included, but ONLY if the corresponding CS files / run IDs already
exist and are passed through the grid.

Recommended use
---------------
Start with the evaluation-only mode (outer CV, no nested selection). This gives an
honest out-of-fold estimate for each fixed configuration.

If later you want to *select one best configuration*, enable nested CV.

@Author: L-F-S
"""

from __future__ import annotations

import itertools
import math
import os
import sys
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, f1_score, precision_score, recall_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2] if len(HERE.parents) >= 3 else HERE.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from conf import CS_DIR, DATA_DIR, DISEASE,IC50_DRUG_COLLAPSE_METHOD,IC50_EFF_TH,\
    IC50_ONLY, IMG_DIR, LINCS_METADATA_PATH,LOGS_DIR, cell_line,cell_line_run_name,\
    chembl_val_log_filename,cs_filename,cs_log_filename,cs_mith,cs_on_LM,disease_run_name,\
    ic50_file, selected_cs_run_id
from logger import append_run_metadata 
from loader import load_IC50, load_drug_rankings
from final_step_correlations_IC50_CS import add_pert_iname_to_cs, classify_ic50_vs_cs,\
    collapse_profiles_to_drug, plot_binchen_fig3_style, resolve_cs_run_id

# -----------------------------------------------------------------------------
# config object
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class EvalConfig:
    """
    Container for one hyperparameter configuration.
    """
    name: str
    disease: str
    cell_line: str
    cs_out_dir: Path
    cs_run_id: str
    cs_filename: str

    cs_drug_colname: str = "pert_iname"
    ic50_drug_colname: str = "pert_iname"
    cs_score_col: str = "connectivity_score"
    ic50_score_col: str = "standard_value"
    cs_drug_collapse_method: str = "best"
    ic50_drug_collapse_method: str = "median"
    ic50_only: bool = True

    cs_threshold: float = -1.5
    ic50_threshold: float = 10.0

    lincs_metadata_path: Path = LINCS_METADATA_PATH
    ic50_file: Path = ic50_file


# -----------------------------------------------------------------------------
# utilities
# -----------------------------------------------------------------------------


def expand_grid(param_grid):
    """
    Create cartesian product of a parameter grid.

    param_grid: dict
        {param_name: list_of_values}

    returns: list of dicts
        Each dict is one combination of parameters
    """
    keys = list(param_grid.keys())
    values = [param_grid[k] for k in keys]
    out = []
    for combo in itertools.product(*values):
        out.append(dict(zip(keys, combo)))
    return out



def pick_n_splits(y, max_splits=5, min_splits=2):
    """
    Choose valid number of CV splits based on class counts.
    """
    y = np.asarray(y)
    _, counts = np.unique(y, return_counts=True)
    n_splits = int(min(max_splits, counts.min()))
    if n_splits < min_splits:
        raise ValueError("Too few samples for stratified CV")
    return n_splits



def resolve_cs_file(config):
    """
    Resolve CS file path from config.
    """
    if config.cs_filename:
        return Path(config.cs_out_dir) / config.cs_filename
    return Path(config.cs_out_dir) / f"{config.cs_run_id}.tsv"


# -----------------------------------------------------------------------------
# data preparation
# -----------------------------------------------------------------------------


def load_merged_for_config(config):
    """
    Load + collapse + merge CS and IC50 for a given config.
    """
    cs_file = resolve_cs_file(config)

    dr = load_drug_rankings(config.cs_out_dir, filename=cs_file.name)
    dr = add_pert_iname_to_cs(dr, config.lincs_metadata_path)

    dr = collapse_profiles_to_drug(
        dr,
        score_col=config.cs_score_col,
        drug_col=config.cs_drug_colname,
        how=config.cs_drug_collapse_method,
    )

    ic50, _ = load_IC50(
        config.ic50_file,
        config.disease,
        cell_line=config.cell_line,
        IC50_ONLY=config.ic50_only,
    )

    ic50 = collapse_profiles_to_drug(
        ic50,
        score_col=config.ic50_score_col,
        drug_col=config.ic50_drug_colname,
        how=config.ic50_drug_collapse_method,
    )

    merged = pd.merge(
        dr,
        ic50,
        left_on=config.cs_drug_colname,
        right_on=config.ic50_drug_colname,
        how="inner",
    )

    merged["actual_positive"] = (
        pd.to_numeric(merged[config.ic50_score_col], errors="coerce")
        <= config.ic50_threshold
    ).astype(int)

    merged = merged.dropna(subset=[config.cs_score_col, config.ic50_score_col])
    merged["drug_name"] = merged[config.cs_drug_colname].astype(str)

    return merged.reset_index(drop=True)


# -----------------------------------------------------------------------------
# CV evaluation
# -----------------------------------------------------------------------------


def compute_fold_metrics(test_df, config):
    """
    Compute classification + correlation metrics on one test fold.
    """
    classified_df, pr = classify_ic50_vs_cs(
        test_df,
        score_col=config.cs_score_col,
        ic50_col=config.ic50_score_col,
        cs_threshold=config.cs_threshold,
        ic50_threshold=config.ic50_threshold,
    )

    y_true = classified_df["actual_positive"].astype(int)
    y_pred = classified_df["predicted_positive"].astype(int)

    f1 = f1_score(y_true, y_pred, zero_division=0)

    return {
        "precision": pr["precision"],
        "recall": pr["recall"],
        "f1": f1,
    }, classified_df



def run_cv_for_config(config, n_splits, n_repeats=10, random_state=42):
    """
    Run repeated stratified CV for one config.
    """
    merged = load_merged_for_config(config)
    y = merged["actual_positive"].values

    splitter = RepeatedStratifiedKFold(
        n_splits=n_splits,
        n_repeats=n_repeats,
        random_state=random_state,
    )

    oof = []
    folds = []

    for i, (_, test_idx) in enumerate(splitter.split(merged, y)):
        test_df = merged.iloc[test_idx]
        metrics, classified = compute_fold_metrics(test_df, config)

        metrics["fold"] = i
        folds.append(metrics)

        classified["fold"] = i
        oof.append(classified)

    oof = pd.concat(oof, ignore_index=True)
    folds = pd.DataFrame(folds)

    return oof, folds


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------


def make_default_config_grid():
    """
    Define hyperparameter grid.
    """
    run_id = resolve_cs_run_id(
        cs_runs_tsv=cs_log_filename,
        disease_run_id=disease_run_name,
        drug_run_id=cell_line_run_name,
        cs_on_LM=cs_on_LM,
        mith=cs_mith,
        selected_cs_run_id=selected_cs_run_id,
    )

    grid = {
        "cs_drug_collapse_method": ["best", "median"],
        "ic50_drug_collapse_method": ["median"],
        "cs_threshold": [-1.8, -1.6, -1.5],
        "ic50_only": [True, False],
        "cell_line": [cell_line, None],
    }

    configs = []
    for i, params in enumerate(expand_grid(grid)):
        name = f"cfg_{i}"
        configs.append(
            EvalConfig(
                name=name,
                disease=DISEASE,
                cell_line=params["cell_line"],
                cs_out_dir=Path(CS_DIR),
                cs_run_id=run_id,
                cs_filename=f"{run_id}.tsv",
                cs_drug_collapse_method=params["cs_drug_collapse_method"],
                ic50_drug_collapse_method=params["ic50_drug_collapse_method"],
                cs_threshold=params["cs_threshold"],
                ic50_only=params["ic50_only"],
            )
        )

    return configs

#%%
if __name__ == "__main__":
    configs = make_default_config_grid()

    first_df = load_merged_for_config(configs[0])
    n_splits = pick_n_splits(first_df["actual_positive"].values)

    for cfg in configs:
        print("Running:", cfg.name)
        oof, folds = run_cv_for_config(cfg, n_splits)
        print(folds.mean())
