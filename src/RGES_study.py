# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:00:20 2025

@author: L-F-S
"""

import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import time
from conf import IMG_DIR
from connectivity_score import montecarlo_connectivity, calculate_RGES
import matplotlib.pyplot as plt
import  matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#%% 
# PARAMETERS:
# connectivity score function:
score_type='bin_chen'#'lamb'#'sirota'#
# score_type='evil_twin'#'lamb'#'sirota'#

# Length of drug gene signature:
r=1000#14812#9892

# Length of disease up regulated signature:
# set to sup=int(r/2) to sample drug spectrum
s_up=int(r/2)#10#7884#4455

# Length of disease up regulated signature:
# set to sup=int(r/2) to sample drug spectrum
s_down=int(r/2)#10#6928#3954

# Number of random iterations:
n_iterations=10000

#%% sample connectivity score distribution over of random rankings of r genes:
    
start=time.time()
css=np.round(montecarlo_connectivity(s_up, s_down, r, n_iterations, score_type=score_type), 3)
print('iterations:', n_iterations,'\nr=',r,'\ntime:', time.time()-start)

#%%

def nice_hist(x, title,  imgname,xlabel='', ylabel='', save=False):
    # plt.rcParams.update({'font.size': 22})
    plt.figure()
    plt.hist(x, bins=100, color='black')
    plt.grid(linestyle='--')
    plt.title(title)
    plt.ylabel(ylabel, fontsize=18)
    plt.xlabel(xlabel, fontsize=18)
    if save:
        plt.savefig(imgname+'.pdf')
    return
    
# imgname='up'+str(s_up)+'down'+str(s_down)+'r'§+str(r)+score_type
imgname='hist_nice'
xlabel='punteggio di connettività'
ylabel='numerosità'
nice_hist(css, '', imgname, xlabel,ylabel)
#%%
