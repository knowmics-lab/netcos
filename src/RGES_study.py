# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:00:20 2025

@author: L-F-S
"""

import os

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
r=10000#14812#9892

# Length of disease up regulated signature:
# set to s_up=int(r/2) to sample drug spectrum
s_up=int(r/2)#10#7884#4455

# Length of disease down regulated signature:
# set to s_down=int(r/2) to sample drug spectrum
s_down=int(r/2)#10#6928#3954

# Number of random iterations:
n_iterations=10000000

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
#%% Split css into n random drug samplings, and pllot correlations between them
# and FC based drug ranking and netocos based drug ranking

# Todo: maybe bring onto notebook or make different notebook about it? 
import numpy as np
from scipy import stats
from conf import CS_OUT
import mataplotlib.pyplot as plt

n_random_drug_rankings=1000 #  n_iterations must be multiple of this number
random_drug_rankings = np.reshape(np.array(css),(n_random_drug_rankings,int(n_iterations/n_random_drug_rankings)))  # c

#%% build correlations
def calc_correlations(measured_ranking, random_rankings):
    KSs=[]
    for random_rank in random_rankings:
        # spearmans.append(stats.spearmanr(measured_ranking,random_rank))
        KSs.append(stats.ks_2samp(measured_ranking,random_rank))
    return  zip(*KSs)

mith_cs_data=pd.read_csv(CS_OUT+'mith_connectivity_score.tsv', sep='\t')
NetCos_KSs, NetCos_KSs_pval =  calc_correlations(mith_cs_data.connectivity_score, random_drug_rankings)
DEG_cs_data=pd.read_csv(CS_OUT+'DEG_connectivity_score.tsv', sep='\t')
DEG_KSs, DEG_KSs_pval =  calc_correlations(DEG_cs_data.connectivity_score, random_drug_rankings)

#%% 
plt.boxplot(x=[NetCos_KSs, DEG_KSs], vert=True, tick_labels = ['PF','FC'])
plt.title('KS test distribution between calculated CS and 1000 random sampling of CS function')
plt.savefig(IMG_DIR+'KS_test_Boxplot_vs_random_samplings')
# da questi risultati sembra peggio sembra che siamo piu casuali noi
#
#%% Altra idea, invece di vedere quanto il mnostro drug ranking si distanzi
# da un ranking di disease casuale, quindi confrontando le nostre 3k drug FC oppure
# drug PF contro tipo 10k ranking casuali di geni. Quindi prendere le FC vere
# e confrontarle contro drug rankings, a questo punto potrei fare delle correlazioni
# xke cosi ho un drug ranking sulle drug esistenti vere e proprie vs un disease c
# casuale. Questa sicuramente da mettere in un altro notebook

