#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:34:00 2020

@author: kate
"""
#Package imports
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#GAMES imports
import Settings
from Solvers import calcChi2, calcRsq, solveSingle
from Run import addNoise, solveAll
from Saving import createFolder

#Define settings
conditions_dictionary, initial_params_dictionary, data_dictionary, HBS_info = Settings.init()
run = conditions_dictionary["run"]
full_path = conditions_dictionary["directory"]
model = conditions_dictionary["model"]
data = conditions_dictionary["data"]
error = data_dictionary["error"]
exp_data = data_dictionary["exp_data"]
location = conditions_dictionary["location"]

#Set style file
if location == 'desktop':
    plt.style.use('C://Users/Katie_Dreyer/Google_Drive/Documents/Leonard_Lab/HBS_Modeling/HBS_GAMES/paper.mplstyle.py')
    
elif location == 'laptop':
    plt.style.use('/Users/kdreyer/Desktop/Github/HBS_GAMES/paper.mplstyle.py')

elif location == 'quest':
    plt.style.use('/home/ksd844/HBS_GAMES/paper.mplstyle.py')


# def saveRefData(data_):
    
#     ''' 
#     Purpose: Save reference training data
        
#     Inputs: 
#         data: list of floats, each float describes reporter expression for a different condition 
#             (length = # datapoints )
           
#     Output:
#         df_ref: dataframe containing the reference training data in a structure compatible with 
#             downstream GAMES code 
        
#     Files: 
#         REFERENCE TRAINING DATA.xlsx (dataframe containing reference training data)
    
#     '''
#     df_ref = pd.DataFrame(data_, columns = ['DsRE2'])
    

        
#     #Save results
#     filename = 'REFERENCE TRAINING DATA.xlsx'

#     with pd.ExcelWriter(filename) as writer:  # doctest: +SKIP
#         df_ref.to_excel(writer, sheet_name = 'Ref')
    
#     return df_ref


def generateRefData(p_ref):
    
    ''' 
    Purpose: Generate reference training data (with and without noise) based on parameter set, p_ref
        
    Inputs: 
        p_ref: list of floats, each float corresponds to a reference parameter value 
                (length = # free parameters). Parameter labels defined in init() 
                (in Settings.py) and SolveAll() (in Run.py).
           
    Output: None
    
    Figures: Model_Fit.svg (plot showing sim vs. training data)
    '''
    
    os.chdir(full_path)
    name = ('All_States_' + model)
    sub_folder_name = './' + name
    createFolder('./' + sub_folder_name)
    os.chdir('./' + sub_folder_name)
    
    #Solve for simulation data
    t_hox, solutions = solveAll(p_ref, exp_data, 'plotting', model, 'all states')
    
    exp_1a = exp_data[:6]
    exp_4b = exp_data[6:13]
    exp_4c = exp_data[13:]
    
    err_1a = error[:6]
    err_4b = error[6:13]
    err_4c = error[13:]
    
    t_1a = [0, 24, 48, 72, 96]
    t_4b = [0, 24, 48, 72, 96, 120]
    t_4c = t_4b 
    
    fig = plt.figure(figsize = (6.5,2.5))
    fig.subplots_adjust(wspace=0.2)
    ax1 = plt.subplot(131)   
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133)
    
    sky_blue = [i/255 for i in [86, 180, 233]]
    colors = [sky_blue, 'gray']
    marker_ = 'o'
    linestyles = ['dotted', 'dotted', 'dotted']

    
    ax1.errorbar(t_1a, exp_1a[:-1], color = colors[0], marker = marker_, yerr = err_1a[:-1], 
                 fillstyle = 'none', linestyle = 'none',capsize = 2, label = '1% O2 training data')
    ax1.errorbar(t_1a[0], exp_1a[-1], color = colors[1], marker = marker_, yerr = err_1a[-1], 
                 fillstyle = 'none', linestyle = 'none',capsize = 2, label = '21% O2 training data')
    
    ax2.errorbar(t_4b, exp_4b[:-1], color = colors[0], marker = marker_, yerr = err_4b[:-1], 
                 fillstyle = 'none', linestyle = 'none',capsize = 2, label = '1% O2 training data')
    ax2.errorbar(t_4b[0], exp_4b[-1], color = colors[1], marker = marker_, yerr = err_4b[-1], 
                 fillstyle = 'none', linestyle = 'none',capsize = 2, label = '21% O2 training data')
    
    ax3.errorbar(t_4c, exp_4c[:-1], color = colors[0], marker = marker_, yerr = err_4c[:-1], 
                 fillstyle = 'none', linestyle = 'none',capsize = 2, label = '1% O2 training data')
    ax3.errorbar(t_4c[0], exp_4c[-1], color = colors[1], marker = marker_, yerr = err_4c[-1], 
             fillstyle = 'none', linestyle = 'none',capsize = 2, label = '21% O2 training data')
  
     #Plot simulated data for the best case parameter set
    ax1.plot(t_hox[:26], solutions[0][:-1], marker = None, label = '1% O2 model fit', 
             linestyle = linestyles[0], color = colors[0])
    ax1.plot(t_hox[0], solutions[0][-1], marker = None, label = '21% O2 model fit', 
             linestyle = linestyles[0], color = colors[1])
    
    ax2.plot(t_hox, solutions[1][:-1], marker = None, label = '1% O2 model fit', 
             linestyle = linestyles[0], color = colors[0])
    ax2.plot(t_hox[0], solutions[1][-1], marker = None, label = '21% O2 model fit', 
             linestyle = linestyles[0], color = colors[1])
    
    ax3.plot(t_hox, solutions[2][:-1], marker = None, label = '1% O2 model fit', 
             linestyle = linestyles[0], color = colors[0])
    ax3.plot(t_hox[0], solutions[2][-1], marker = None, label = '21% O2 model fit', 
             linestyle = linestyles[0], color = colors[1])
    
    #Set x and y labels and ylim
    ax1.set_xlabel('Time Post-Plating (hours)')
    ax1.set_ylabel('Relative DsRed Expression')
    ax1.set_ylim(ax2.get_ylim())
    ax1.set_title('Simple HBS')
    ax1.set_box_aspect(1)
    ax1.legend()
    
    ax2.set_xlabel('Time Post-Plating (hours)')
    ax2.set_ylabel('Relative DsRed Expression')
    ax2.set_title('HIF1a Feedback HBS')
    ax2.set_box_aspect(1)
    ax2.legend()
    
    ax3.set_xlabel('Time Post-Plating (hours)')
    ax3.set_ylabel('Relative DsRed Expression')
    ax3.set_ylim(ax2.get_ylim())
    ax3.set_title('HIF2a Feedback HBS')
    ax3.set_box_aspect(1)
    ax3.legend()

    plt.savefig('./Model_Fit.svg', bbox_inches="tight")
    
    return
            
    # #Calculate chi2 between reference training data with and without noise
    # chi2 = calcChi2(noise_solutions, norm_solutions, error)
    # print('chi2(p_ref) = ' + str(round(chi2, 4)))
    
    # #Save reference dataframe as an excel sheet
    # saveRefData(noise_solutions)

p_ref = [1.0, 1.0, 1.0e-1, 1.0, 1.0, 1.0e-1, 1.0, 1.0, 1.0e-1, 1.0e-2,
         1.0e-2, 1.0, 1.0e-2, 1.0e-2]

generateRefData(p_ref)

    
