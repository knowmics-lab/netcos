# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:41:29 2025

@author: L-F-S

Run MITHrIL batch on LINCS drug signatures and log the run.
"""

import subprocess
import time
import socket
from pathlib import Path
from datetime import datetime

import pandas as pd

from conf import (
    MITH_APP, MITH_OUT_DRUG, MITH_IN_DRUG, mith_threads, mith_organism,
    mith_input_file, LOGS_DIR, cell_line, cell_line_run_name, DISEASE, disease_run_name,
    landmark_disease, landmark_drug, LM_flag_disease, LM_flag_drug,
)


def append_mithril_run_metadata(metadata_path, row_dict):
    """
    Append one row to a MITHrIL run log TSV, creating the file if needed.
    """
    metadata_path = Path(metadata_path)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    row_df = pd.DataFrame([row_dict])

    if metadata_path.exists():
        old_df = pd.read_csv(metadata_path, sep="\t")
        out_df = pd.concat([old_df, row_df], ignore_index=True)
    else:
        out_df = row_df

    out_df.to_csv(metadata_path, sep="\t", index=False)


def run_mithril_batch(
    mith_input_file,
    MITH_APP,
    MITH_IN_DRUG,
    MITH_OUT_DRUG,
    organism="hsa",
    n_thread="30",
    verbose=True,
    printc=False,
    log_run=True,
):
    """
    Run MITHrIL batch on a LINCS drug-signature batch input file.
    """

    in_file_path = Path(MITH_IN_DRUG) / mith_input_file
    out_dir = Path(MITH_OUT_DRUG)

    command = [
        "java", "-jar", str(MITH_APP),
        "mithril-batch", "-reactome", "-p", "-customize-pathway-matrix",
        "-seed", "1234",
        "-inversion-factory", "fast-cpu",
        "-multiplication-factory", "fast-cpu",
        "-t", str(n_thread),
        "-organism", organism,
        "-i", str(in_file_path),
        "-o", str(out_dir),
        "-p",
    ]

    if verbose:
        command.append("-verbose")

    if printc:
        print(" ".join(command))
        return None

    out_dir.mkdir(parents=True, exist_ok=True)

    n_files_before = len(list(out_dir.glob("*")))

    start = time.time()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    stdout_lines = []
    stderr_lines = []

    for line in process.stdout:
        print(line, end="")
        stdout_lines.append(line)
    for line in process.stderr:
        print(line, end="")
        stderr_lines.append(line)

    process.wait()
    elapsed_sec = time.time() - start

    n_files_after = len(list(out_dir.glob("*")))
    n_new_files = n_files_after - n_files_before
    n_perturbation_files = len(list(out_dir.glob("*.perturbations.txt")))
    n_output_txt_files = len([p for p in out_dir.glob("*.txt") if not p.name.endswith(".perturbations.txt")])

    if log_run:
        timestamp = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
        mithril_run_id = f"{timestamp}__{cell_line_run_name}"

        metadata_row = {
            "mithril_run_id": mithril_run_id,
            "timestamp": timestamp,
            "hostname": socket.gethostname(),
            "step_name": "step6_drug_run_MITHrIL",
            "run_scope": "drug_batch",
            "drug_run_id": cell_line_run_name,
            "disease_run_id": disease_run_name,
            "disease_symbol": DISEASE,
            "cell_line": cell_line,
            "filter_disease_signature_by_landmark_genes_pre_mith": int(bool(landmark_disease)),
            "filter_drug_signature_by_landmark_genes_pre_mith": int(bool(landmark_drug)),
            "LM_flag_disease": LM_flag_disease,
            "LM_flag_drug": LM_flag_drug,
            "mith_input_file": str(in_file_path),
            "mith_output_dir": str(out_dir),
            "mith_input_filename": str(mith_input_file),
            "mith_app": str(MITH_APP),
            "mith_organism": organism,
            "mith_threads": int(n_thread),
            "verbose": int(bool(verbose)),
            "return_code": process.returncode,
            "elapsed_sec": elapsed_sec,
            "n_output_files_before": n_files_before,
            "n_output_files_after": n_files_after,
            "n_new_output_files": n_new_files,
            "n_perturbation_files_after": n_perturbation_files,
            "n_main_txt_files_after": n_output_txt_files,
            "stdout_n_lines": len(stdout_lines),
            "stderr_n_lines": len(stderr_lines),
        }

        append_mithril_run_metadata(
            Path(LOGS_DIR) / "mithril_drug_runs.tsv",
            metadata_row
        )

    return process.returncode


if __name__ == "__main__":
    from sys import argv

    in_name = mith_input_file
    if len(argv) > 1:
        in_name = argv[1]

    run_mithril_batch(
        in_name,
        MITH_APP,
        MITH_IN_DRUG,
        MITH_OUT_DRUG,
        n_thread=mith_threads,
        organism=mith_organism,
        printc=False,
    )
