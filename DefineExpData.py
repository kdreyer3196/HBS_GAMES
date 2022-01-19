#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kateldray
"""
#Package imports
import numpy as np

def defineExp(data, df_ref, df_err):
    ''' 
    Purpose: Import the experimental data from an (.xlsx sheet) and 
        structure the data to be compatible with downstream LMFIT code.
    
    Inputs:
        data: a string defining the data identity
            (= 'ligand dose response only' or 'ligand dose response and DBD dose response')
        df_ref: a dataframe containing the reference data
      
        
    Outputs:
         x: a list of lists defining the component doses for each datapoint. 
             Each list has 3 items - [DBD dose, AD dose, ligand dose] - and there 
             are the same number of lists as number of datapoints.
         data: a list of floats defining the normalized reporter expression values 
             for each datapoint (length = # datapoints)
         error: list of floats defining the normalized reporter expression error 
             values for each datapoint (length = # datapoints)
        
  '''

    exp_lists = df_ref.values.tolist()
    err_lists = df_err.values.tolist()

    
    if data == 'hypox only':
        exp1a = np.append(exp_lists[0][:5], exp_lists[1][3])
        exp4b = np.append(exp_lists[2], exp_lists[3][3])
        exp4c = np.append(exp_lists[4], exp_lists[5][3])
        
        err1a = np.append(err_lists[0][:5], err_lists[1][3])
        err4b = np.append(err_lists[2], err_lists[3][3])
        err4c = np.append(err_lists[4], err_lists[5][3])
        
        
    elif data == 'all':
        
        exp1a = np.concatenate((exp_lists[0][:5], exp_lists[1][3:]))
        exp4b = np.concatenate((exp_lists[2], exp_lists[3][3:]))
        exp4c = np.concatenate((exp_lists[4], exp_lists[5][3:]))

        err1a = np.concatenate((err_lists[0][:5], err_lists[1][3:]))
        err4b = np.concatenate((err_lists[2], err_lists[3][3:]))
        err4c = np.concatenate((err_lists[4], err_lists[5][3:]))
    
    
    data = np.concatenate((exp1a, exp4b, exp4c))
    error = np.concatenate((err1a, err4b, err4c))

       
    return data, error
