# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:41:29 2025

@author: L-F-S
"""

import subprocess
from conf import MITH_APP, MITH_OUT_DRUG, MITH_IN_DRUG

def run_mithril_batch( MITH_APP, MITH_IN_DRUG, MITH_OUT_DRUG, organism='hsa',
                 n_thread="30", verbose=True, printc=True):
    '''
    MITHrIL batch parameters:
    -organism hsa_v2023_03 organism version used in this experiment. For most 
            up-to-date version, use 'hsa'. for a list of all available
            organisms, type:
            java -jar MITH_APP organisms

    -reactome: also use reactome network
    -p a flag (unlike non-batch version) to create perturbation file
    -customize-pathway-matrix : calcola alcone cose come mithril 2. va messo
    -inversion-factory fast-cpu : accelereate calculation with fast-cpu library
    -multiplication-factory fast-cpu : accelereate calculation with fast-cpu library
    -t : number of threads
    -o : output directory (unlike non-batch version, where is output name)'''
    
    command = [
        "java", "-jar", MITH_APP,
        "mithril-batch", "-reactome", "-p", "-customize-pathway-matrix",
        "-seed", "1234", "-inversion-factory", "fast-cpu", "-multiplication-factory", "fast-cpu",
        "-t", n_thread, "-organism", organism,
        "-i", f"{MITH_IN_DRUG}LINCS_metanalysis.mi ",
        "-o", f"{MITH_OUT_DRUG}",
        "-p"
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

run_mithril_batch(MITH_APP, MITH_IN_DRUG, MITH_OUT_DRUG, n_thread="30")
