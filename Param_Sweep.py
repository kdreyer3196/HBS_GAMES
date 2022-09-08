# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 18:42:23 2022

@author: Katie_Dreyer
"""
# =============================================================================
# IMPORTS
# =============================================================================
#Import external packages/functions
import datetime
import os
import multiprocessing as mp
import math
from math import log10
import warnings
import pandas as pd
import numpy as np
from lmfit import Model, Parameters
from openpyxl import load_workbook
import matplotlib.pyplot as plt 
import seaborn as sns
import cycler
from math import sqrt

#Import GAMES functions
from Solvers import *
import Settings
from Run import solveAll, solvePar 

warnings.filterwarnings("ignore")

#Unpack conditions from Settings.py
conditions_dictionary, initial_params_dictionary, data_dictionary, HBS_info = Settings.init()
run = conditions_dictionary["run"]
model = conditions_dictionary["model"] 
data = conditions_dictionary["data"]
num_cores = conditions_dictionary["num_cores"]
n_initial_guesses = conditions_dictionary["n_initial_guesses"]
problem_free = conditions_dictionary["problem"]
full_path = conditions_dictionary["directory"]
n_search = conditions_dictionary["n_search"]
n_search_parameter_estimation= conditions_dictionary["n_search"]
n_search_pem_eval = conditions_dictionary["n_search_pem_eval"]
modules = conditions_dictionary["modules"] 
confidence_interval = conditions_dictionary["confidence_interval"] 
fit_params = problem_free['names']
bounds = problem_free['bounds']
num_vars = problem_free['num_vars']
p_all = conditions_dictionary["p_all"] 
p_ref = conditions_dictionary["p_ref"] 
p_labels_free = conditions_dictionary["p_labels_free"] 
p_ref_free = conditions_dictionary["p_ref_free"] 
real_param_labels_all = conditions_dictionary["real_param_labels_all"] 
num_datasets_pem_eval = conditions_dictionary["num_datasets_pem_eval"] 
param_index_PL = conditions_dictionary["param_index_PL"]
problem_all_params = conditions_dictionary["problem_all_params"]
all_param_labels = problem_all_params['names']
param_labels = list(initial_params_dictionary.keys())
init_params = list(initial_params_dictionary.values())
real_param_labels_free = conditions_dictionary["real_param_labels_free"]
data_type = data_dictionary["data_type"]
exp_data = data_dictionary["exp_data"]
exp_data_original = data_dictionary["exp_data"]
error = data_dictionary["error"]
location = conditions_dictionary["location"]

save_internal_states_flag = True

#Set style file
if location == 'desktop':
    plt.style.use('C://Users/Katie_Dreyer/Desktop/Github/HBS_GAMES/paper.mplstyle.py')
    parallelization = 'no' 
        
elif location == 'laptop':
    plt.style.use('/Users/kdreyer/Desktop/Github/HBS_GAMES/paper.mplstyle.py')
    parallelization = 'no' 

elif location == 'quest':
    plt.style.use('/home/ksd844/HBS_GAMES/paper.mplstyle.py')
    parallelization = 'yes' 

#Define pO2 conditions
O2_range = [7.6, 138]

def get_param_sweeps(p_ref, orders_of_mag):
    
    param_sweeps = {}
    
    for index, param in enumerate(p_ref):
        param_label = real_param_labels_all[index]
        
        param_sweeps[param_label] = []
        
        
        for exp in orders_of_mag:
            param_new = param*10**exp
            
            param_sweeps[param_label].append(param_new)  
        
    return param_sweeps

def get_param_lists(p_ref, param_sweeps):
    
    param_lists = {}
    
    for index, param in enumerate(p_ref):
        
        param_label = real_param_labels_all[index]
        
        param_lists[param_label + ' sweep'] = []
            
        for p_new in param_sweeps[param_label]:
            
            new_params = p_ref.copy()
    
            new_params[index] = p_new
                        
            param_lists[param_label + ' sweep'].append(new_params)
            
    return param_lists

# def run_param_sims(param_lists):
    
#     for param_label, param_sets in param_lists.items():

#         for params in param_sets:
#             t_hox, SS_hox_1a, SS_hox_4b, SS_hox_4c, norm = solveAll(params, exp_data, 'plotting', model, ' ')

        

# def plot_param_sweep():
    
    
# os.chdir(full_path)

orders_of_mag = [-3, -2, -1, 0, 1, 2, 3]

param_sweeps = get_param_sweeps(p_ref, orders_of_mag)
param_lists = get_param_lists(p_ref, param_sweeps)
print(param_lists)