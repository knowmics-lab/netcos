# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:41:29 2025

@author: L-F-S
"""

import subprocess
from conf import DISEASE, MITH_APP, MITH_OUT_DISEASE, MITH_IN_DISEASE

def run_mithril(DISEASE, MITH_APP, MITH_IN_DISEASE, MITH_OUT_DISEASE,
                organism='hsa_v2023_03', n_thread="30", verbose=True, printc=True): 
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
    
    
    command = [
        "java", "-jar", MITH_APP,
        "mithril", "-reactome", "-p", "-customize-pathway-matrix",
        "-seed", "1234", "-inversion-factory", "fast-cpu", "-multiplication-factory", "fast-cpu",
        "-t", n_thread, "-i", f"{MITH_IN_DISEASE}{DISEASE}_signature_gene_id.mi",
        "-organism", organism,
        "-o", f"{MITH_OUT_DISEASE}{DISEASE}_mith3.output.txt",
        "-p", f"{MITH_OUT_DISEASE}{DISEASE}_mith3.perturbations.txt"
    ]

    if verbose:
        command.append("-verbose")
    if printc:
        print(' '.join(command))
        return 
    else:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        for line in process.stdout:
            print(line, end="")
        for line in process.stderr:
            print(line, end="")
        
        process.wait()
    
        return process.returncode

run_mithril(DISEASE, MITH_APP, MITH_IN_DISEASE, MITH_OUT_DISEASE)