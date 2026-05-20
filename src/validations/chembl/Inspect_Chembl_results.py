# -*- coding: utf-8 -*-
"""
Inspect top ChEMBL validation runs.
Creates a compact TSV summary of the top runs per cell line.
"""

import os
import sys

HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

import pandas as pd

from conf import chembl_val_log_filename, diseases_of, LOGS_DIR


##############################################################################
# PARAMETERS
##############################################################################

SORT_BY = "spearman_r"
TOP_N = 10
ASCENDING = False

OUT_FILE = LOGS_DIR / f"top_{TOP_N}_chembl_validation_runs_by_{SORT_BY}.tsv"

COLUMNS_TO_SHOW = [
    "cell_line",
    "disease_run_id",
    "spearman_r",
    "spearman_pval",
    "precision",
    "recall",
    "n_unique_drugs",
    "landmark_drug",
    "landmark_disease",
    "CS_METHOD",
    "cs_mith",
    "cs_on_LM",
    "CS_ON_PATHWAYS",
    "LINCS_drug_collapse_method",
    "IC50_drug_collapse_method",
    "ic50_only",
    "cell_line_IC50",
    "cs_threshold",
    "ic50_threshold",
    "cs_run_id",
]


##############################################################################
# LOAD
##############################################################################

df = pd.read_csv(chembl_val_log_filename, sep="\t")

if SORT_BY not in df.columns:
    raise ValueError(
        f"Column '{SORT_BY}' not found.\n"
        f"Available columns:\n{list(df.columns)}"
    )


##############################################################################
# COLLECT TOP RUNS
##############################################################################

out_rows = []

for cell_line, disease in diseases_of.items():

    cell_df = df[df["cell_line_IC50"] == cell_line].copy()

    if cell_df.empty:
        print(f"No rows found for {cell_line}")
        continue

    top_df = (
        cell_df
        .sort_values(SORT_BY, ascending=ASCENDING)
        .head(TOP_N)
        .reset_index(drop=True)
    )
    
    top_df["rank"] = range(1, len(top_df) + 1)
    top_df["disease"] = disease
    top_df["cell_line"] = cell_line
    
    out_rows.append(top_df)


summary_df = pd.concat(out_rows, ignore_index=True)

existing_cols = [col for col in COLUMNS_TO_SHOW if col in summary_df.columns]
summary_df = summary_df[existing_cols]

summary_df.to_csv(OUT_FILE, sep="\t", index=False)


##############################################################################
# PRINT COMPACT SUMMARY
##############################################################################

print("\nTop ChEMBL validation runs")
print("=" * 80)
PRINT_COLUMNS = [
    "disease_run_id",
    "spearman_r",
    "precision",
    "recall",
    "n_unique_drugs",
    "CS_METHOD",
    "cs_mith",
    "cs_on_LM",
    "CS_ON_PATHWAYS",
    "LINCS_drug_collapse_method",
    "IC50_drug_collapse_method",
]

for cell_line in diseases_of.keys():

    print(f"\n{cell_line}")
    print("-" * 80)

    small = summary_df[summary_df["cell_line"] == cell_line]
    
    print(small[PRINT_COLUMNS].to_string(index=False))

print(f"\nSaved TSV:\n{OUT_FILE}")
print("\nDone.")