#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert MITHrIL pathway summary output (*.output.txt) into NetCos pathway-level CS input.

This is the pathway-level analogue of step7:
- for disease:  MITHrIL output -> CS input
- for drugs:    MITHrIL outputs -> CS input(s)

It uses MITHrIL's pathway summary table (*.output.txt), NOT the perturbation table.
The pathway effect size used for CS is the signed 'Corrected Accumulator'.

Output columns
--------------
- pathway_id
- pathway_name
- perturbation
- p.value
"""

import os
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd

HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO_ROOT)

from conf import  DISEASE, MITH_OUT_DISEASE, MITH_OUT_DRUG, CS_IN_DISEASE, CS_IN_DRUG, LOGS_DIR
from logger import append_run_metadata


STEP_NAME = "step2b_7b_mith_pathway_to_cs_input"
TIMESTAMP = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
LOG_FILE = Path(LOGS_DIR) / f"{DISEASE}_step2b_7b_runs.tsv"


REQUIRED_COLS = ["# Pathway Id",  "Pathway Name",  "Corrected Accumulator", "Adjusted pValue"]


def load_mithril_pathway_output(path):
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")
    return df


def build_pathway_cs_input(df):
    out = df.loc[:, REQUIRED_COLS].copy()
    out = out.rename(columns={"# Pathway Id": "ID"})

    out["ID"] = out["ID"].astype(str)
    out["Pathway Name"] = out["Pathway Name"].astype(str)
    out["Corrected Accumulator"] = pd.to_numeric(out["Corrected Accumulator"], errors="coerce")
    out["Adjusted pValue"] = pd.to_numeric(out["Adjusted pValue"], errors="coerce")

    out = out.dropna(subset=["ID", "Corrected Accumulator", "Adjusted pValue"]).copy()
    out = out.drop_duplicates(subset=["ID"], keep="first")
    out = out.sort_values("ID").reset_index(drop=True)

    return out


def convert_one_file(in_file, out_file):
    df = load_mithril_pathway_output(in_file)
    out = build_pathway_cs_input(df)

    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_file, sep="\t", index=False)

    print(f"[OK] {in_file} -> {out_file} ({len(out)} pathways)")

    return {"input_file": str(in_file), "output_file": str(out_file), "n_input_rows": df.shape[0],"n_output_rows": out.shape[0],  }


def build_disease_paths():
    in_file = os.path.join(MITH_OUT_DISEASE, f"{DISEASE}_mith3.output.txt")
    out_file = os.path.join(CS_IN_DISEASE, f"{DISEASE}_mith_pathways.tsv")
    return in_file, out_file


def build_drug_paths():
    mith_dir = Path(MITH_OUT_DRUG)
    cs_dir = Path(CS_IN_DRUG)
    cs_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for in_file in sorted(mith_dir.glob("*.output.txt")):
        stem = in_file.name.replace(".output.txt", "")
        out_file = cs_dir / f"{stem}_pathways.tsv"
        jobs.append((str(in_file), str(out_file)))
    return jobs


# ---------------------------
# RUNNERS
# ---------------------------

def run_disease():
    in_file, out_file = build_disease_paths()

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"Disease MITHrIL output not found: {in_file}")

    meta = convert_one_file(in_file, out_file)

    log_row = {
        "timestamp": TIMESTAMP,
        "step": STEP_NAME,
        "mode": "disease",
        "disease": DISEASE,
        **meta,
    }

    append_run_metadata(LOG_FILE, log_row)


def run_drugs():
    jobs = build_drug_paths()

    if not jobs:
        raise FileNotFoundError(f"No drug MITHrIL outputs found in: {MITH_OUT_DRUG}")

    for in_file, out_file in jobs:
        meta = convert_one_file(in_file, out_file)

        log_row = {
            "timestamp": TIMESTAMP,
            "step": STEP_NAME,
            "mode": "drug",
            "disease": DISEASE,
            **meta,
        }

        append_run_metadata(LOG_FILE, log_row)
#%%

if __name__ == "__main__":
    run_disease()
    run_drugs()