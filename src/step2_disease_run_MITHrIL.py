# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:41:29 2025

@author: L-F-S
"""

import subprocess
from conf import DISEASE, MITH_APP, MITH_OUT_DISEASE, MITH_IN_DISEASE

def run_mithril(DISEASE, MITH_APP, MITH_IN_DISEASE, MITH_OUT_DISEASE, verbose=True):
    command = [
        "java", "-jar", MITH_APP,
        "mithril", "-reactome", "-p", "-customize-pathway-matrix",
        "-seed", "1234", "-inversion-factory", "fast-cpu", "-multiplication-factory", "fast-cpu",
        "-t", "30", "-i", f"{MITH_IN_DISEASE}{DISEASE}_signature_gene_id.mi",
        "-o", f"{MITH_OUT_DISEASE}{DISEASE}_mith3.output.txt",
        "-p", f"{MITH_OUT_DISEASE}{DISEASE}_mith3.perturbations.txt"
    ]
    if verbose:
        command.append("-verbose")
    
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    for line in process.stdout:
        print(line, end="")
    for line in process.stderr:
        print(line, end="")
    
    process.wait()

    return process.returncode

run_mithril(DISEASE, MITH_APP, MITH_IN_DISEASE, MITH_OUT_DISEASE)
