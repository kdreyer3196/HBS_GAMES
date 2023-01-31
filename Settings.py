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
    folder_name = '230126_ModelD2_Modules1-2'
    
    # =============================================================================
    # 2. Define conditions dictionary
    # =============================================================================
    #Initialize conditions dictionary
    #Items that you might want to change
    conditions_dictionary = {}
    conditions_dictionary["model"] = 'model_D' #'model_A', 'model_B', 'model_C', or 'model_D'
    conditions_dictionary["modules"] = [1,2] #[1,2,3] or [1,2] or [2,3] or [1] or [2] or [3] or [] for test only
    conditions_dictionary["n_search"] = 1000
    conditions_dictionary["n_initial_guesses"] = 100
    conditions_dictionary["confidence_interval"] = .99 
    conditions_dictionary["num_cores"] = 28
    conditions_dictionary["num_datasets_pem_eval"] = 3
    conditions_dictionary["n_search_pem_eval"] = 1000
    conditions_dictionary["param_index_PL"] = 'all' #'all' or index of p (int)
    conditions_dictionary["data"] = 'hypox only'
    conditions_dictionary["location"] = 'quest' #'desktop', 'laptop', 'quest'
    
    # =============================================================================
    # 3. Define free parameters and bounds
    # =============================================================================
    #Set list of all potentially free parameters
    
    p_ref = [36, 1.0, 1.0e-2, 0.1, 0.1, 1.0, 1.0e-1, 1.0e-2, 1.0e-2, 1.0]
                 
    [t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH] = p_ref
    
    if conditions_dictionary["model"] == 'model_A':
        t_HAF = 0.0
        k_txn2 = 1.0
        k_dHAF = 0.0
        k_bHS = 0.0

        p_ref = [
            t_HAF, 
            k_txn2, 
            k_dHAF, 
            k_bHS, 
            k_bHH, 
            k_txnH, 
            k_dH1R, 
            k_dH1P, 
            k_dHP, 
            k_txnBH
        ]
        
        real_param_labels_free = [
            'k_bHH',
            'k_txnH',
            'k_dH1R',
            'k_dH1P', 
            'k_dHP', 
            'k_txnBH'
        ]
    
    if conditions_dictionary["model"] == 'model_B':
        # t_HAF = 0.0
        # k_txn2 = 1.0
        # k_dHAF = 0.0
        k_bHS = 0.0

        p_ref = [
            t_HAF, 
            k_txn2, 
            k_dHAF, 
            k_bHS, 
            k_bHH, 
            k_txnH, 
            k_dH1R, 
            k_dH1P, 
            k_dHP, 
            k_txnBH
        ]
        
        real_param_labels_free = [
            't_HAF',
            'k_txn2',
            'k_dHAF',
            'k_bHH',
            'k_txnH',
            'k_dH1R',
            'k_dH1P', 
            'k_dHP', 
            'k_txnBH'
        ]

    if conditions_dictionary["model"] == 'model_C':
        t_HAF = 0.0
        k_txn2 = 1.0
        k_dHAF = 0.0
        # k_bHS = 0.0

        p_ref = [
            t_HAF, 
            k_txn2, 
            k_dHAF, 
            k_bHS, 
            k_bHH, 
            k_txnH, 
            k_dH1R, 
            k_dH1P, 
            k_dHP, 
            k_txnBH
        ]
        
        real_param_labels_free = [
            'k_bHS',
            'k_bHH',
            'k_txnH',
            'k_dH1R',
            'k_dH1P', 
            'k_dHP', 
            'k_txnBH'
        ]

    if conditions_dictionary["model"] == 'model_D':
        # t_HAF = 0.0
        # k_txn2 = 1.0
        # k_dHAF = 0.0
        # k_bHS = 0.0
        # k_bHH = 0.0

        p_ref = [
            t_HAF, 
            k_txn2, 
            k_dHAF, 
            k_bHS, 
            k_bHH, 
            k_txnH, 
            k_dH1R, 
            k_dH1P, 
            k_dHP, 
            k_txnBH
        ]
        
        real_param_labels_free = [
            't_HAF',
            'k_txn2',
            'k_dHAF',
            'k_bHS',
            'k_bHH',
            'k_txnH',
            'k_dH1R',
            'k_dH1P', 
            'k_dHP', 
            'k_txnBH'
        ]

    p_all = p_ref

    #Define parameter labels (real and general)
    
    #real labels for p_ref and p_all     
    
    real_param_labels_all = [
        't_HAF',
        'k_txn2', 
        'k_dHAF',
        'k_bHS', 
        'k_bHH',
        'k_txnH', 
        'k_dH1R',
        'k_dH1P',
        'k_dHP',
        'k_txnBH'
    ]

    #general labels for p_ref and p_all

    p_labels_all = [
        'p1',
        'p2', 
        'p3', 
        'p4', 
        'p5', 
        'p6', 
        'p7', 
        'p8', 
        'p9',
        'p10'
        ]
    
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

    if real_param_labels_free[0] == 't_HAF':
        bounds_log[0] = [1.3802, 1.8573]
            
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
    
    #Model A
    HBS_info['HBS_1aA'] = {}
    HBS_info['HBS_1aA']['# states'] = 10
    HBS_info['HBS_1aA']['state names'] = ['HAFR', 'HAFP', 'aHIF',
                                          'HIF1R', 'HIF1P', 'HIF2R',
                                          'HIF2P', 'HIF2P*', 'DSRE2R',
                                          'DSRE2P']
    
    HBS_info['HBS_4bA'] = {}
    HBS_info['HBS_4bA']['# states'] = 10
    HBS_info['HBS_4bA']['state names'] = ['HAFR', 'HAFP', 'aHIF',
                                          'HIF1R', 'HIF1P', 'HIF2R',
                                          'HIF2P', 'HIF2P*', 'DSRE2R',
                                          'DSRE2P']

    HBS_info['HBS_4cA'] = {}
    HBS_info['HBS_4cA']['# states'] = 10
    HBS_info['HBS_4cA']['state names'] = ['HAFR', 'HAFP', 'aHIF',
                                          'HIF1R', 'HIF1P', 'HIF2R',
                                          'HIF2P', 'HIF2P*', 'DSRE2R',
                                          'DSRE2P']

    #Model B
    HBS_info['HBS_1aB'] = {}
    HBS_info['HBS_1aB']['# states'] = 10
    HBS_info['HBS_1aB']['state names'] = ['HAFR', 'HAFP', 'aHIF',
                                          'HIF1R', 'HIF1P', 'HIF2R',
                                          'HIF2P', 'HIF2P*', 'DSRE2R',
                                          'DSRE2P']
    
    HBS_info['HBS_4bB'] = {}
    HBS_info['HBS_4bB']['# states'] = 10
    HBS_info['HBS_4bB']['state names'] = ['HAFR', 'HAFP', 'aHIF',
                                          'HIF1R', 'HIF1P', 'HIF2R',
                                          'HIF2P', 'HIF2P*', 'DSRE2R',
                                          'DSRE2P']

    HBS_info['HBS_4cB'] = {}
    HBS_info['HBS_4cB']['# states'] = 10
    HBS_info['HBS_4cB']['state names'] = ['HAFR', 'HAFP', 'aHIF',
                                          'HIF1R', 'HIF1P', 'HIF2R',
                                          'HIF2P', 'HIF2P*', 'DSRE2R',
                                          'DSRE2P']

    #Model C
    HBS_info['HBS_1aC'] = {}
    HBS_info['HBS_1aC']['# states'] = 13
    HBS_info['HBS_1aC']['state names'] = ['HAFR', 'HAFP', 'SUMOR', 'SUMOP',
                                           'HAFS', 'aHIF', 'HIF1R', 'HIF1P',
                                           'HIF2R', 'HIF2P', 'HIF2P*', 'DSRE2R',
                                           'DSRE2P']
    
    HBS_info['HBS_4bC'] = {}
    HBS_info['HBS_4bC']['# states'] = 13
    HBS_info['HBS_4bC']['state names'] = ['HAFR', 'HAFP', 'SUMOR', 'SUMOP',
                                           'HAFS', 'aHIF', 'HIF1R', 'HIF1P',
                                           'HIF2R', 'HIF2P', 'HIF2P*', 'DSRE2R',
                                           'DSRE2P']

    HBS_info['HBS_4cC'] = {}
    HBS_info['HBS_4cC']['# states'] = 13
    HBS_info['HBS_4cC']['state names'] = ['HAFR', 'HAFP', 'SUMOR', 'SUMOP',
                                           'HAFS', 'aHIF', 'HIF1R', 'HIF1P',
                                           'HIF2R', 'HIF2P', 'HIF2P*', 'DSRE2R',
                                           'DSRE2P']

    #Model D
    HBS_info['HBS_1aD'] = {}
    HBS_info['HBS_1aD']['# states'] = 13
    HBS_info['HBS_1aD']['state names'] = ['HAFR', 'HAFP', 'SUMOR', 'SUMOP',
                                           'HAFS', 'aHIF', 'HIF1R', 'HIF1P',
                                           'HIF2R', 'HIF2P', 'HIF2P*', 'DSRE2R',
                                           'DSRE2P']
    
    HBS_info['HBS_4bD'] = {}
    HBS_info['HBS_4bD']['# states'] = 13
    HBS_info['HBS_4bD']['state names'] = ['HAFR', 'HAFP', 'SUMOR', 'SUMOP',
                                           'HAFS', 'aHIF', 'HIF1R', 'HIF1P',
                                           'HIF2R', 'HIF2P', 'HIF2P*', 'DSRE2R',
                                           'DSRE2P']

    HBS_info['HBS_4cD'] = {}
    HBS_info['HBS_4cD']['# states'] = 13
    HBS_info['HBS_4cD']['state names'] = ['HAFR', 'HAFP', 'SUMOR', 'SUMOP',
                                           'HAFS', 'aHIF', 'HIF1R', 'HIF1P',
                                           'HIF2R', 'HIF2P', 'HIF2P*', 'DSRE2R',
                                           'DSRE2P']

    return conditions_dictionary, initial_params_dictionary, data_dictionary, HBS_info
