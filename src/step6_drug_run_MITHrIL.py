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
import pandas as pd
from conf import cell_line, MITH_APP, MITH_OUT_DRUG, MITH_IN_DRUG, mith_input_file, \
    mith_batch_threads, mith_organism, LOGS_DIR, cell_line_run_name, landmark_drug
from logger import append_run_metadata

#%% logger functions

def get_first_output_file(output_dir):
    output_dir = Path(output_dir)

    perturb_files = sorted(output_dir.glob("*.perturbations.txt"))

    if not perturb_files:
        return None

    return perturb_files[0]

def count_mithril_drug_outputs(output_dir):
    """
    Count drug-level perturbation output files produced by mithril-batch.
    Assumes one *.perturbations.txt file per drug/signature output.
    """
    output_dir = Path(output_dir)
    return len(list(output_dir.glob("*.perturbations.txt")))

def get_n_unique_output_ids(first_output_file):
    """
    Count unique entities in the MITHrIL perturbation output, using the 'Gene Id' column.
    This includes genes, miRNAs, metabolites, etc.
    """
    first_output_file = Path(first_output_file)

    df = pd.read_csv(first_output_file, sep="\t")

    gene_id_col = "Gene Id"
    if gene_id_col not in df.columns:
        raise ValueError(f"Column '{gene_id_col}' not found in {first_output_file}")

    ids = df[gene_id_col].dropna().astype(str).str.strip()
    ids = ids[ids != ""]

    return ids.nunique()

def get_n_lines_from_first_output(output_dir):
    """
    Return the number of gene rows in the first perturbation output file.
    Assumes MITHrIL output is a tab-separated text file with one header row.
    Returns None if no perturbation files are present.
    """
    output_dir = Path(output_dir)
    perturb_files = sorted(output_dir.glob("*.perturbations.txt"))

    if len(perturb_files) == 0:
        return None, None

    first_file = perturb_files[0]

    # count lines minus header
    with open(first_file, "r", encoding="utf-8", errors="replace") as f:
        n_lines = sum(1 for _ in f)

    n_genes = max(n_lines - 1, 0)
    return n_genes, str(first_file)


#%%
def run_mithril_batch(mith_input_file, MITH_APP, MITH_IN_DRUG, MITH_OUT_DRUG, organism='hsa',
                 n_thread="30", verbose=True, printc=True, run_mith=True):
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
    
    return_code=123
    if run_mith:
        start = time.time()
        process = subprocess.Popen(command)
    
        process.wait()
        return_code = process.returncode
        elapsed_sec = time.time() - start

    # n_files_after = len(list(out_dir.glob("*")))
    # n_new_files = n_files_after - n_files_before
    # n_perturbation_files = len(list(out_dir.glob("*.perturbations.txt")))
    # n_output_txt_files = len([p for p in out_dir.glob("*.txt") if not p.name.endswith(".perturbations.txt")])

    n_files_after = len(list(out_dir.glob("*")))
    n_new_files = n_files_after - n_files_before
    
    n_drugs_output = count_mithril_drug_outputs(out_dir)
    n_output_txt_files = len([p for p in out_dir.glob("*.txt") if not p.name.endswith(".perturbations.txt")])
    
    n_lines_first_output, first_output_file = get_n_lines_from_first_output(out_dir)
    n_gene_ids=get_n_unique_output_ids(first_output_file)

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
        "return_code": return_code,
        "elapsed_sec": elapsed_sec,
        "n_output_files_before": n_files_before,
        "n_output_files_after": n_files_after,
        "n_new_output_files": n_new_files,
        # "n_perturbation_files_after": n_perturbation_files,
        "n_main_txt_files_after": n_output_txt_files,
        
        "n_drugs_output": n_drugs_output,
        "n_genes_first_output": n_gene_ids,
        "first_output_file_used_for_gene_count": first_output_file,
    }
    
    append_run_metadata(Path(LOGS_DIR) / "mithril_drug_runs.tsv", metadata_row)
    return 

#%% 
if __name__ == "__main__":
    from sys import argv
    if len(argv)>1:
        mith_input_file = argv[1]
    run_mithril_batch(mith_input_file, MITH_APP, MITH_IN_DRUG, MITH_OUT_DRUG, \
                      n_thread=mith_batch_threads, organism=mith_organism, printc=False, run_mith=False)
