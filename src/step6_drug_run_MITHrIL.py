# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:41:29 2025

@author: L-F-S
"""

import subprocess
import time
import socket
from pathlib import Path
from datetime import datetime
from conf import cell_line, MITH_APP, MITH_OUT_DRUG, MITH_IN_DRUG, mith_input_file\
    mith_batch_threads, mith_organism, LOGS_DIR, cell_line_run_name, landmark_drug
from logger import append_run_metadata

def run_mithril_batch(mith_input_file, MITH_APP, MITH_IN_DRUG, MITH_OUT_DRUG, organism='hsa_2025',
                 n_thread="30", verbose=True, printc=True):
    '''
    MITHrIL batch parameters:
    -organism: hsa_v2023_03 organism version used in this experiment. For most 
            up-to-date version, use default 'hsa'. for a list of all available
            organisms, type:
            java -jar MITH_APP organisms

    -reactome: also use reactome network
    -seed: random seed for replicability
    -p a flag (unlike non-batch version) to create perturbation file
    -customize-pathway-matrix : calcola alcone cose come mithril 2. va messo
    -inversion-factory fast-cpu : accelereate calculation with fast-cpu library
    -multiplication-factory fast-cpu : accelereate calculation with fast-cpu library
    -t : number of threads
    -o : output directory (unlike non-batch version, where is output name)'''
    
    in_file_path = MITH_IN_DRUG / mith_input_file
    out_dir = Path(MITH_OUT_DRUG)

    command = [
        "java", "-jar", MITH_APP,
        "mithril-batch", "-reactome", "-p", "-customize-pathway-matrix",
        "-seed", "1234", "-inversion-factory", "fast-cpu", "-multiplication-factory", "fast-cpu",
        "-t", str(n_thread), "-organism", organism,
        "-i", f"{in_file_path}",
        "-o", str(out_dir),
        "-p"
    ]

    if verbose:
        command.append("-verbose")
    
    if printc:
        print(' '.join(command))
        return 

    # esle: run MITHrIL        
    out_dir.mkdir(parents=True, exist_ok=True)
    n_files_before = len(list(out_dir.glob("*")))

    start = time.time()
    process = subprocess.Popen(command)

    process.wait()
    elapsed_sec = time.time() - start

    n_files_after = len(list(out_dir.glob("*")))
    n_new_files = n_files_after - n_files_before
    n_perturbation_files = len(list(out_dir.glob("*.perturbations.txt")))
    n_output_txt_files = len([p for p in out_dir.glob("*.txt") if not p.name.endswith(".perturbations.txt")])

    # Log run
    timestamp = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    mithril_run_id = f"{timestamp}__{cell_line_run_name}"

    metadata_row = {
        "mithril_run_id": mithril_run_id,
        "timestamp": timestamp,
        "hostname": socket.gethostname(),
        "step_name": "step6_drug_run_MITHrIL",
        "run_scope": "drug_batch",
        "drug_run_id": cell_line_run_name,
        "cell_line": cell_line,
        "filter_drug_signature_by_landmark_genes_pre_mith": int(bool(landmark_drug)),
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
    }
    
    append_run_metadata(Path(LOGS_DIR) / "mithril_drug_runs.tsv", metadata_row)
    return process.returncode

#%% 
if __name__ == "__main__":
    from sys import argv
    if len(argv)>1:
        mith_input_file = argv[1]
    run_mithril_batch(mith_input_file, MITH_APP, MITH_IN_DRUG, MITH_OUT_DRUG, \
                      n_thread=mith_batch_threads, organism=mith_organism, printc=True)
