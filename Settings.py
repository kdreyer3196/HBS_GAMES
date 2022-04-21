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
    folder_name = '220420_Model1_PEM_new_bounds'
    
    # =============================================================================
    # 2. Define free parameters and bounds
    # =============================================================================
    #Set list of all potentially free parameters

    # Model 0
    # p_ref = [7.91, 1.43e-4, 1.51e-2, 3.91e-3, 1.08e-2, 1.07, 1.15e-2, 9.96e-2,
             # 0.967] 
    
    # Model 1
    p_ref = [1.43e-4, 1.51e-2, 1.08e-2, 1.15e-2, 1.0e-2, 1.0e-2]
    
    p_all = p_ref

    #Define parameter labels (real and general)
    
    # Model 0
    # real_param_labels_all = ['k_prod', 'k_pMARS', 'k_out', 'k_dMARS',
                             # 'k_dreg1', 'k_tln2', 'k_dreg2', 'k_dreg3',
                             # 'deg_ratio']
                             
    # Model 1
    real_param_labels_all = ['k_pMARS', 'k_dMARS', 'k_dreg1', 'k_dO2',
                             'k_dreg2', 'k_act']
    
    #real labels for p_ref and p_all
    
    # Model 0
    # p_labels_all = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9'] #general param labels
    
    # Model 1
    p_labels_all = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6']
    
    #if a param in real_param_labels_all is not included in realParamLabels_free,
    #it is fixed at the value set in p_all
    
    # Model 0
    # real_param_labels_free = ['k_prod', 'k_pMARS', 'k_out', 'k_dMARS',
                             # 'k_dreg1', 'k_tln2', 'k_dreg2', 'k_dreg3',
                             # 'deg_ratio']
    
    # Model 1
    real_param_labels_free = ['k_pMARS', 'k_dMARS', 'k_dreg1', 'k_dO2',
                             'k_dreg2', 'k_act']  #real labels for free params

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
    
    # Model 0
    # bounds_kprod = [-4, 1]
    # bounds_kpMARS = [-4, 1]
    # bounds_kout = [-4, 2]
    # bounds_kdMARS = [-4, 1]
    # bounds_kdreg1 = [-4, 1]
    # bounds_ktln2 = [-1, 2]
    # bounds_kdreg2 = [-4, 1]
    # bounds_kdreg3 = [-4, 1]
    # bounds_degratio = [-1, 1]
    
    # bounds_log = [bounds_kprod, bounds_kpMARS, bounds_kout, bounds_kdMARS, 
    #               bounds_kdreg1, bounds_ktln2, bounds_kdreg2, bounds_kdreg3, 
    #               bounds_degratio]
    
    # Model 1
    bounds_kpMARS = [-4, 2]
    bounds_kdMARS = [-4, 2]
    bounds_kdreg1 = [-5, 2]
    bounds_kdO2 = [-4, 2]
    bounds_kdreg2 = [-5, 2]
    bounds_kact = [-4, 2]
    
    bounds_log = [bounds_kpMARS, bounds_kdMARS, bounds_kdreg1, bounds_kdO2, 
                  bounds_kdreg2, bounds_kact]
        
    #Define the parameter estimation problem (free parameters for this run only)
    problem = {'num_vars': num_free_params,  #set free parameters and bounds
               'names': p_labels_free, 
               'bounds': bounds_log} #bounds are in log scale
    # =============================================================================
    # 3. Define conditions dictionary
    # =============================================================================
    #Initialize conditions dictionary
    #Items that you might want to change
    conditions_dictionary = {}
    conditions_dictionary["model"] = 'model 1' #'model 0', 'model 1'
    conditions_dictionary["modules"] = [2] #[1,2,3] or [1,2] or [2,3] or [1] or [2] or [3] or [] for test only
    conditions_dictionary["n_search"] = 20
    conditions_dictionary["n_initial_guesses"] = 5
    conditions_dictionary["confidence_interval"] = .99 
    conditions_dictionary["num_cores"] = 14
    conditions_dictionary["num_datasets_pem_eval"] = 3
    conditions_dictionary["n_search_pem_eval"] = 20
    conditions_dictionary["param_index_PL"] = 'all' #'all' or index of p (int)
    conditions_dictionary["data"] = 'hypox only'

    #Define parameter labels
    #Items that you likely will not change
    full_path = makeMainDir(folder_name)
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
    
    #laptop file path
    # path_exp = ('/Users/kdreyer/Google Drive/My Drive/Documents/Leonard_Lab/HBS_Modeling/' +
            # 'Resources and Notes/Experimental_Data/20210512_Exp2_Analysis_for_Katie.xlsx')
    
    #desktop file path
    # path_exp = ('C://Users/Katie_Dreyer/Google_Drive/Documents/Leonard_Lab/HBS_Modeling/' +
            # 'Resources and Notes/Experimental_Data/20210512_Exp2_Analysis_for_Katie.xlsx')
    
    #QUEST file path
    path_exp = ('/home/ksd844/HBS_GAMES/Exp_Data.xlsx')
    
    df_ref = pd.read_excel(path_exp, sheet_name='Exp_Data_Norm', header=0, index_col=[0,1], engine='openpyxl')
    df_err = pd.read_excel(path_exp, sheet_name='Exp_Error_Norm', header=0, index_col=[0,1], engine='openpyxl')
    exp_data, error = defineExp(conditions_dictionary["data"], df_ref, df_err)
    data_dictionary["exp_data"] = exp_data
    data_dictionary["error"] = error
    data_dictionary["data_type"] = ''
    
    HBS_info = {}
    
    # MODEL 0
    HBS_info['HBS_1a'] = {}
    HBS_info['HBS_1a']['# states'] = 8
    HBS_info['HBS_1a']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4b'] = {}
    HBS_info['HBS_4b']['# states'] = 9
    HBS_info['HBS_4b']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'sHIF1R', 'HIF1P',
                                         'HIF2R', 'HIF2P', 'DSRed2R', 'DSRed2P']
    
    HBS_info['HBS_4c'] = {}
    HBS_info['HBS_4c']['# states'] = 9
    HBS_info['HBS_4c']['state names'] = ['MARS0', 'MARS', 'HIF1R', 'HIF1P', 'HIF2R',
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
    
    return conditions_dictionary, initial_params_dictionary, data_dictionary, HBS_info
