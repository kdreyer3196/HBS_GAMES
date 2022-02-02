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

#local file path
#plt.style.use('C://Users/Katie_Dreyer/Google_Drive/Documents/Leonard_Lab/HBS_Modeling/HBS_GAMES/paper.mplstyle.py')

#QUEST file path
plt.style.use('/home/ksd844/HBS_GAMES/paper.mplstyle.py')


def saveRefData(data_):
    
    ''' 
    Purpose: Save reference training data
        
    Inputs: 
        data: list of floats, each float describes reporter expression for a different condition 
            (length = # datapoints )
           
    Output:
        df_ref: dataframe containing the reference training data in a structure compatible with 
            downstream GAMES code 
        
    Files: 
        REFERENCE TRAINING DATA.xlsx (dataframe containing reference training data)
    
    '''
    df_ref = pd.DataFrame(data_, columns = ['DsRE2'])
    

        
    #Save results
    filename = 'REFERENCE TRAINING DATA.xlsx'

    with pd.ExcelWriter(filename) as writer:  # doctest: +SKIP
        df_ref.to_excel(writer, sheet_name = 'Ref')
    
    return df_ref


def generateRefData(p_ref):
    
    ''' 
    Purpose: Generate reference training data (with and without noise) based on parameter set, p_ref
        
    Inputs: 
        p_ref: list of floats, each float corresponds to a reference parameter value 
                (length = # free parameters). Parameter labels defined in init() 
                (in Settings.py) and SolveAll() (in Run.py).
           
    Output: None
    
    Figures: TRAINING DATA.svg (plot showing reference training data)
    '''
    
    os.chdir(full_path)
    sub_folder_name = './REFERENCE DATA'
    createFolder('./' + sub_folder_name)
    os.chdir('./' + sub_folder_name)
    
    #Solve for simulation data
    norm_solutions, chi2 = solveAll(p_ref, exp_data)
    
    ref_1a = norm_solutions[:6]
    ref_4b = norm_solutions[6:13]
    ref_4c = norm_solutions[13:]
    
    #Add technical error
    noise_solutions = addNoise(norm_solutions, 0)

    noise_ref_1a = noise_solutions[:6]
    noise_ref_4b = noise_solutions[6:13]
    noise_ref_4c = noise_solutions[13:]
        
    t_1a = [0, 24, 48, 72, 96]
    t_4b = [0, 24, 48, 72, 96, 120]
    t_4c = t_4b 

    #Plot 
    fig = plt.figure(figsize = (9,3))
    fig.subplots_adjust(wspace=0.2)
    ax1 = plt.subplot(131)   
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133)
    
    colors = ['dimgrey', 'black']
    linestyles = ['dotted', 'dotted', 'dotted']
    
      #Plot simulated data for the best case parameter set
    ax1.plot(t_1a, ref_1a[:-1], marker = None, label = '1% O2 true', 
             linestyle = linestyles[0], color = colors[0])
    ax1.plot(t_1a[0], ref_1a[-1], marker = None, label = '21% O2 true', 
             linestyle = linestyles[0], color = colors[1])
    
    ax2.plot(t_4b, ref_4b[:-1], marker = None, label = '1% O2 true', 
             linestyle = linestyles[0], color = colors[0])
    ax2.plot(t_4b[0], ref_4b[-1], marker = None, label = '21% O2 true', 
             linestyle = linestyles[0], color = colors[1])
    
    ax3.plot(t_4c, ref_4c[:-1], marker = None, label = '1% O2 true', 
             linestyle = linestyles[0], color = colors[0])
    ax3.plot(t_4c[0], ref_4c[-1], marker = None, label = '21% O2 true', 
             linestyle = linestyles[0], color = colors[1])
    ######################################################################
    ax1.plot(t_1a, noise_ref_1a[:-1], marker = 'o', label = '1% O2 noise', 
             linestyle = 'none', markerSize = 6, fillstyle = 'none', color = colors[0])
    ax1.plot(t_1a[0], noise_ref_1a[-1], marker = 'o', label = '21% O2 noise', 
             linestyle = 'none', markerSize = 6, fillstyle = 'none', color = colors[1])
    
    ax2.plot(t_4b, noise_ref_4b[:-1], marker = 'o', label = '1% O2 noise', 
             linestyle = 'none', markerSize = 6, fillstyle = 'none', color = colors[0])
    ax2.plot(t_4b[0], noise_ref_4b[-1], marker = 'o', label = '21% O2 noise', 
             linestyle = 'none', markerSize = 6, fillstyle = 'none', color = colors[1])
    
    ax3.plot(t_4c, noise_ref_4c[:-1], marker = 'o', label = '1% O2 noise', 
             linestyle = 'none', markerSize = 6, fillstyle = 'none', color = colors[0])
    ax3.plot(t_4c[0], noise_ref_4c[-1], marker = 'o', label = '21% O2 noise', 
             linestyle = 'none', markerSize = 6, fillstyle = 'none', color = colors[1])
    
    #Set x and y labels
    ax1.set_xlabel('Time Post-Plating (hours)')
    ax1.set_ylabel('Relative DsRE2 Expression')
    ax1.set_title('Simple HBS')
    ax1.legend()
    
    ax2.set_xlabel('Time Post-Plating (hours)')
    ax2.set_ylabel('Relative DsRE2 Expression')
    ax2.set_title('HIF1a Feedback HBS')
    ax2.legend()
    
    ax3.set_xlabel('Time Post-Plating (hours)')
    ax3.set_ylabel('Relative DsRE2 Expression')
    ax3.set_title('HIF2a Feedback HBS')
    ax3.legend()
    
    plt.savefig('./REFERENCE TRAINING DATA.svg', dpi = 600)
        
    #Calculate chi2 between reference training data with and without noise
    chi2 = calcChi2(noise_solutions, norm_solutions, error)
    print('chi2(p_ref) = ' + str(round(chi2, 4)))
    
    #Save reference dataframe as an excel sheet
    saveRefData(noise_solutions)

p_ref = [7.91, 1.43e-4, 1.51e-2, 3.91e-3, 1.08e-2, 1.07, 1.15e-2, 9.96e-2, 0.967]    
generateRefData(p_ref)



    
