#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:12:12 2020

@author: kate
"""
#Package imports
from math import log10
import pandas as pd
import os

#GAMES imports
from Saving import makeMainDir
from DefineExpData import defineExp


def init():
    '''
    Purpose: Define conditions for simulations, called by other functions to import conditions
    
    Input: None
    
    Outputs: 3 dictionaries - conditions_dictionary, initial_params_dictionary, data_dictionary
    
    '''
    # =============================================================================
    # 1. Define and create folder for saving results
    # =============================================================================
    #This will be the name of the run-specific results folder. 
    folder_name = '220427_Model2A_Modules1-2'
    
    # =============================================================================
    # 2. Define conditions dictionary
    # =============================================================================
    #Initialize conditions dictionary
    #Items that you might want to change
    conditions_dictionary = {}
    conditions_dictionary["model"] = 'model 2A' #'model 0', 'model 0A', 'model 0B', 'model 1', 'model 1A', 'model 2', 'model 2A', 'model 3'
    conditions_dictionary["modules"] = [1] #[1,2,3] or [1,2] or [2,3] or [1] or [2] or [3] or [] for test only
    conditions_dictionary["n_search"] = 1000
    conditions_dictionary["n_initial_guesses"] = 100
    conditions_dictionary["confidence_interval"] = .99 
    conditions_dictionary["num_cores"] = 14
    conditions_dictionary["num_datasets_pem_eval"] = 3
    conditions_dictionary["n_search_pem_eval"] = 1000
    conditions_dictionary["param_index_PL"] = 'all' #'all' or index of p (int)
    conditions_dictionary["data"] = 'hypox only'
    conditions_dictionary["location"] = 'quest' #'desktop', 'laptop', 'quest'
    
    # =============================================================================
    # 3. Define free parameters and bounds
    # =============================================================================
    #Set list of all potentially free parameters
    
    p_ref = [7.91, 1.51e-2, 1.43e-4, 3.91e-3, 1.08e-2, 1.15e-2, 9.96e-2,
             1.15e-2, 1.0e-2, 1.07, 0.967]
    
    [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio] = p_ref
    
    if conditions_dictionary["model"] == 'model 0':
        k_dHP = 0
        k_aH2P = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_pM1', 'k_dM1', 'k_dH1R',
                                  'k_dH1P', 'k_pH2R', 'k_tln2', 'deg_ratio']
        
    elif conditions_dictionary["model"] == 'model 0A':
        k_pH2R = 0
        k_dHP = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_pM1', 'k_dM1', 'k_dH1R',
                                  'k_dH1P', 'k_aH2P', 'k_tln2', 'deg_ratio']
        
    elif conditions_dictionary["model"] == 'model 0B':
        k_pH2R = 0
        k_dHP = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_pM1', 'k_dM1', 'k_dH1R',
                                  'k_dH1P', 'k_aH2P', 'k_tln2', 'deg_ratio']

    elif conditions_dictionary["model"] == 'model 1':
        k_pM1 = 0
        k_dM1 = 0
        k_pH2R = 0
        k_tln2 = 0
        deg_ratio = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_dH1R', 'k_dH1P',
                                  'k_dHP', 'k_aH2P']
        
    elif conditions_dictionary["model"] == 'model 1A':
        k_pM1 = 0
        k_dM1 = 0
        k_pH2R = 0
        k_tln2 = 0
        deg_ratio = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_dH1R', 'k_dH1P',
                                  'k_dHP', 'k_aH2P']

    elif conditions_dictionary["model"] == 'model 2':
        k_pH2R = 0
        k_tln2 = 0
        deg_ratio = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_pM1', 'k_dM1', 'k_dH1R',
                                  'k_dH1P', 'k_dHP', 'k_aH2P']
        
    elif conditions_dictionary["model"] == 'model 2A':
        k_pH2R = 0
        k_tln2 = 0
        deg_ratio = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_pM1', 'k_dM1', 'k_dH1R',
                                  'k_dH1P', 'k_dHP', 'k_aH2P']
        
    elif conditions_dictionary["model"] == 'model 3':
        k_pH2R = 0
        deg_ratio = 0
        
        p_ref = [k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP,
                 k_aH2P, k_tln2, deg_ratio]
        
        real_param_labels_free = ['k_pM0', 'k_dM0', 'k_pM1', 'k_dM1', 'k_dH1R',
                                  'k_dH1P', 'k_dHP', 'k_aH2P', 'k_tln2']
    
    p_all = p_ref

    #Define parameter labels (real and general)
    
    #real labels for p_ref and p_all      
    
    real_param_labels_all = ['k_pM0', 'k_dM0', 'k_pM1', 'k_dM1', 'k_dH1R',
                             'k_dH1P', 'k_pH2R', 'k_dHP', 'k_aH2P', 'k_tln2',
                             'deg_ratio']

    #general labels for p_ref and p_all

    p_labels_all = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9',
                    'p10', 'p11']
    
    #if a param in real_param_labels_all is not included in realParamLabels_free,
    #it is fixed at the value set in p_all

    #Change param labels to generalizable param labels
    num_free_params = len(real_param_labels_free)
    initial_params_dictionary = {}
    params = []
    p_ref_free = []
    p_labels_free = []
    for i, value in enumerate(p_all):
        label = real_param_labels_all[i]
        if label in real_param_labels_free:
            initial_params_dictionary[label] = value
            p_labels_free.append(p_labels_all[i])
            params.append(value)
            p_ref_free.append(p_ref[i])
    
    # Set bounds for parameter estimation
    
    bounds_log = []
    for i in range(0, num_free_params):
        min_bound = log10(p_ref_free[i]) - 3
        max_bound = log10(p_ref_free[i]) + 3
        bounds_log.append([min_bound, max_bound])
            
    #Define the parameter estimation problem (free parameters for this run only)
    problem = {'num_vars': num_free_params,  #set free parameters and bounds
               'names': p_labels_free, 
               'bounds': bounds_log} #bounds are in log scale

    
    # =============================================================================
    # 4. Add param info to conditions dictionary
    # =============================================================================
    #Define parameter labels
    #Items that you likely will not change
    full_path = makeMainDir(folder_name, conditions_dictionary["location"])
    conditions_dictionary["real_param_labels_free"] = real_param_labels_free
    problem_all_params = {'num_vars': len(p_labels_all),  
                          'names': p_labels_all, 
                          'bounds': [[]]} 
    conditions_dictionary["real_param_labels_all"] = real_param_labels_all
    conditions_dictionary["p_all"] = p_all
    conditions_dictionary["p_ref_free"] = p_ref_free
    conditions_dictionary["p_labels_free"] = p_labels_free
    conditions_dictionary["p_ref"] = p_ref
    conditions_dictionary["directory"] = full_path
    conditions_dictionary["run"] = '1'
    conditions_dictionary["problem"] = problem
    conditions_dictionary["problem_all_params"] = problem_all_params
    
    # =============================================================================
    # 4. Define data dictionary
    # =============================================================================
    data_dictionary = {}
    
    if conditions_dictionary["location"] == 'desktop':
            path_exp = ('C://Users/Katie_Dreyer/Google_Drive/Documents/Leonard_Lab/HBS_Modeling/' +
            'Resources and Notes/Experimental_Data/20210512_Exp2_Analysis_for_Katie.xlsx')
            
    elif conditions_dictionary["location"] == 'laptop':
            path_exp = ('/Users/kdreyer/Google Drive/My Drive/Documents/Leonard_Lab/HBS_Modeling/' +
            'Resources and Notes/Experimental_Data/20210512_Exp2_Analysis_for_Katie.xlsx')
            
    elif conditions_dictionary["location"] == 'quest':
        path_exp = ('/home/ksd844/HBS_GAMES/Exp_Data.xlsx')
            
    df_ref = pd.read_excel(path_exp, sheet_name='Exp_Data_Norm', header=0, index_col=[0,1], engine='openpyxl')
    df_err = pd.read_excel(path_exp, sheet_name='Exp_Error_Norm', header=0, index_col=[0,1], engine='openpyxl')
    exp_data, error = defineExp(conditions_dictionary["data"], df_ref, df_err)
    data_dictionary["exp_data"] = exp_data
    data_dictionary["error"] = error
    data_dictionary["data_type"] = ''
    
    HBS_info = {}
    
    # MODEL 0
    HBS_info['HBS_1a0'] = {}
    HBS_info['HBS_1a0']['# states'] = 8
    HBS_info['HBS_1a0']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4b0'] = {}
    HBS_info['HBS_4b0']['# states'] = 9
    HBS_info['HBS_4b0']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'sHIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4c0'] = {}
    HBS_info['HBS_4c0']['# states'] = 9
    HBS_info['HBS_4c0']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P', 'HIF2R',
                                         'sHIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    #MODEL 0A
    HBS_info['HBS_1a0A'] = {}
    HBS_info['HBS_1a0A']['# states'] = 8
    HBS_info['HBS_1a0A']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4b0A'] = {}
    HBS_info['HBS_4b0A']['# states'] = 9
    HBS_info['HBS_4b0A']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'sHIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4c0A'] = {}
    HBS_info['HBS_4c0A']['# states'] = 9
    HBS_info['HBS_4c0A']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P', 'HIF2R',
                                         'sHIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    #MODEL 0B
    HBS_info['HBS_1a0B'] = {}
    HBS_info['HBS_1a0B']['# states'] = 8
    HBS_info['HBS_1a0B']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4b0B'] = {}
    HBS_info['HBS_4b0B']['# states'] = 9
    HBS_info['HBS_4b0B']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'sHIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4c0B'] = {}
    HBS_info['HBS_4c0B']['# states'] = 9
    HBS_info['HBS_4c0B']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P', 'HIF2R',
                                         'sHIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    # MODEL 1
    HBS_info['HBS_1a1'] = {}
    HBS_info['HBS_1a1']['# states'] = 7
    HBS_info['HBS_1a1']['state names'] = ['MARS', 'HIF1R', 'HIF1P','HIF2R',
                                          'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4b1'] = {}
    HBS_info['HBS_4b1']['# states'] = 7
    HBS_info['HBS_4b1']['state names'] = ['MARS', 'HIF1R', 'HIF1P','HIF2R',
                                          'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4c1'] = {}
    HBS_info['HBS_4c1']['# states'] = 7
    HBS_info['HBS_4c1']['state names'] = ['MARS', 'HIF1R', 'HIF1P','HIF2R',
                                          'HIF2P', 'DSRed2R', 'DSRed2P']
    
    # MODEL 1A
    HBS_info['HBS_1a1A'] = {}
    HBS_info['HBS_1a1A']['# states'] = 7
    HBS_info['HBS_1a1A']['state names'] = ['MARS', 'HIF1R', 'HIF1P','HIF2R',
                                          'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4b1A'] = {}
    HBS_info['HBS_4b1A']['# states'] = 7
    HBS_info['HBS_4b1A']['state names'] = ['MARS', 'HIF1R', 'HIF1P','HIF2R',
                                          'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4c1A'] = {}
    HBS_info['HBS_4c1A']['# states'] = 7
    HBS_info['HBS_4c1A']['state names'] = ['MARS', 'HIF1R', 'HIF1P','HIF2R',
                                          'HIF2P', 'DSRed2R', 'DSRed2P']
    
    # MODEL 2
    HBS_info['HBS_1a2'] = {}
    HBS_info['HBS_1a2']['# states'] = 8
    HBS_info['HBS_1a2']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    HBS_info['HBS_4b2'] = {}
    HBS_info['HBS_4b2']['# states'] = 8
    HBS_info['HBS_4b2']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    HBS_info['HBS_4c2'] = {}
    HBS_info['HBS_4c2']['# states'] = 8
    HBS_info['HBS_4c2']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    # MODEL 2A
    HBS_info['HBS_1a2A'] = {}
    HBS_info['HBS_1a2A']['# states'] = 8
    HBS_info['HBS_1a2A']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    HBS_info['HBS_4b2A'] = {}
    HBS_info['HBS_4b2A']['# states'] = 8
    HBS_info['HBS_4b2A']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    HBS_info['HBS_4c2A'] = {}
    HBS_info['HBS_4c2A']['# states'] = 8
    HBS_info['HBS_4c2A']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    # MODEL 3
    HBS_info['HBS_1a3'] = {}
    HBS_info['HBS_1a3']['# states'] = 8
    HBS_info['HBS_1a3']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    HBS_info['HBS_4b3'] = {}
    HBS_info['HBS_4b3']['# states'] = 9
    HBS_info['HBS_4b3']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'sHIF1R',
                                          'HIF1P', 'HIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    HBS_info['HBS_4c3'] = {}
    HBS_info['HBS_4c3']['# states'] = 9
    HBS_info['HBS_4c3']['state names'] = ['MARS0', 'MARS1', 'HIF1R', 'HIF1P', 
                                          'HIF2R', 'sHIF2R', 'HIF2P', 'DSRed2R',
                                          'DSRed2P']
    
    return conditions_dictionary, initial_params_dictionary, data_dictionary, HBS_info
