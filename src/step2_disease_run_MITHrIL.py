# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:41:29 2025

@author: L-F-S
"""

import subprocess
from pathlib import Path
import time
from datetime import datetime
import socket
from conf import disease_run_name, MITH_APP, MITH_OUT_DISEASE, MITH_IN_DISEASE,\
    mith_threads, LOGS_DIR, landmark_disease, DISEASE
from logger import append_run_metadata

def run_mithril(disease_run_name, MITH_APP, MITH_IN_DISEASE, MITH_OUT_DISEASE,
                organism='hsa', n_thread="30", verbose=True, printc=True): 
    '''
    MITHrIL parameters:
    -organism hsa_v2023_03 organism version used in this experiment. For most 
            up-to-date version, use 'hsa'. for a list of all available
            organisms, type:
            java -jar MITH_APP organisms

    -reactome: also use reactome network
    -seed: random seed for replicability
    -customize-pathway-matrix : calcola alcone cose come mithril 2. va messo
    -inversion-factory fast-cpu : accelereate calculation with fast-cpu library
    -multiplication-factory fast-cpu : accelereate calculation with fast-cpu library
    -t : number of threads
    -o : full path of output file name
    -p : full path of perturbation output file name
    '''
    
    input_file = Path(MITH_IN_DISEASE) / f"{disease_run_name}_signature_gene_id.mi"
    output_main_file = Path(MITH_OUT_DISEASE) / f"{disease_run_name}_mith3.output.txt"
    output_pert_file = Path(MITH_OUT_DISEASE) / f"{disease_run_name}_mith3.perturbations.txt"

    
    command = [
        "java", "-jar", MITH_APP,
        "mithril", "-reactome", "-p", "-customize-pathway-matrix",
        "-seed", "1234", "-inversion-factory", "fast-cpu", "-multiplication-factory", "fast-cpu",
        "-t", str(n_thread), "-i", str(input_file),
        "-organism", organism,
        "-o", str(output_main_file),
        "-p", str(output_pert_file)
    ]

    if verbose:
        command.append("-verbose")
        
    if printc:
        print(' '.join(command))
        return 
    
    # else: run MITHrIL 
    Path(MITH_OUT_DISEASE).mkdir(parents=True, exist_ok=True)
    
    start = time.time()
    process = subprocess.Popen(command)
    
    process.wait()
    elapsed_sec = time.time() - start
    
    
    main_exists = output_main_file.exists()
    pert_exists = output_pert_file.exists()
    main_size = output_main_file.stat().st_size if main_exists else None
    pert_size = output_pert_file.stat().st_size if pert_exists else None
    
    timestamp = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    mithril_run_id = f"{timestamp}__{disease_run_name}"

    # Log run
    metadata_row = {
        "mithril_run_id": mithril_run_id,
        "timestamp": timestamp,
        "hostname": socket.gethostname(),
        "step_name": "step2_disease_run_MITHrIL",
        "run_scope": "disease",
        "disease_run_id": disease_run_name,
        "disease_symbol": DISEASE,
        "filter_disease_signature_by_landmark_genes_pre_mith": int(bool(landmark_disease)),
        "mith_input_file": str(input_file),
        "mith_output_main_file": str(output_main_file),
        "mith_output_perturbation_file": str(output_pert_file),
        "mith_app": str(MITH_APP),
        "mith_organism": organism,
        "mith_threads": int(n_thread),
        "return_code": process.returncode,
        "elapsed_sec": elapsed_sec,
        "output_main_exists": int(main_exists),
        "output_perturbation_exists": int(pert_exists),
        "output_main_size_bytes": main_size,
        "output_perturbation_size_bytes": pert_size,
    }

    append_run_metadata(Path(LOGS_DIR) / "mithril_disease_runs.tsv", metadata_row)
    
    
    return process.returncode
#%%
if __name__=="__main__":
    from sys import argv
    from conf import disease_run_name
    # insert disease_run_name name 
    # python step2_disease_run_MITHrIL.py disease_run_name
    if len(argv)>1:
        disease_run_name = argv[1]
    run_mithril(disease_run_name, MITH_APP, MITH_IN_DISEASE, MITH_OUT_DISEASE, n_thread=mith_threads,  printc=True)
