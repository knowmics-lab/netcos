#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:54:33 2025

@ author: L-F-S

MITHrIL re-writes all file names that contain characters which are not
alphanumeric, '_', or '-', for safety reasons.
Run this to recuperate original file names, if needed.
Renames mith3 output files by recovering oringinal characters from input files
"""
import os

import pickle
import re
from conf import MITH_IN_DRUG, MITH_OUT_DRUG


#%% check diff between mith3 input and putout filenames
print('mith3 batch output filename bug: some filenames with special characters\
      differ from input to output')

mith_in_drug_names=set([n.split('.')[0] for n in os.listdir(MITH_IN_DRUG)])

def mithril_string_replace(input_string):
    pattern = r'[^a-zA-Z0-9_-]'#r'[^a-zA-Z0-9-]'
    
    # Replace matching characters with a dash
    replaced_string = re.sub(pattern, '-', input_string)
    
    return replaced_string

# print('does not work with this drug, but ok', mithril_string_replace('5\'-guanidinonaltrindole_24h'))

bad_to_good={mithril_string_replace(good):good for good in mith_in_drug_names}

for filename in os.listdir(MITH_OUT_DRUG):
    for bad_name, good_name in bad_to_good.items():
        if bad_name in filename:
            good_filename=filename.replace(bad_name, good_name)
            os.rename(MITH_OUT_DRUG+filename,MITH_OUT_DRUG+good_filename)
            print(filename,'changed to', good_filename, 'after')
            break

#%% 
print('performing final check. Anything printed is not found in output directory.\n should only print LINCS files:\n')
mith_in_drug_names=[n.split('.')[0] for n in os.listdir(MITH_IN_DRUG)]
mith_out_drug_names=[n.split('.')[0] for n in os.listdir(MITH_OUT_DRUG)]

for infile in mith_in_drug_names:
    if not infile in mith_out_drug_names:
        print(infile)