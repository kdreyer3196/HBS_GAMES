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
# from Run import solveAll

#ignore ODEint warnings that clog up the console - user can remove this line if they want to see the warnings
warnings.filterwarnings("ignore")

#Define colors for plotting (colors from a colorblind friendly palette)
black_ = [i/255 for i in [0, 0, 0]]
orange_ = [i/255 for i in [230, 159, 0]]
sky_blue = [i/255 for i in [86, 180, 233]]
pink_ = [i/255 for i in [204, 121, 167]]

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

save_internal_states_flag = False

#Set style file
if location == 'desktop':
    plt.style.use('C://Users/Katie_Dreyer/Desktop/Github/HBS_GAMES/paper.mplstyle.py')
    parallelization = 'no' 
    
elif location == 'laptop':
    plt.style.use('/Users/kdreyer/Documents/Github/HBS_GAMES/paper.mplstyle.py')
    parallelization = 'no' 

elif location == 'quest':
    plt.style.use('/home/ksd844/HBS_GAMES/paper.mplstyle.py')
    parallelization = 'yes' 

def get_df(f_name):
    os.chdir(full_path)
    df = pd.read_excel(f_name, header=0, index_col=0, engine='openpyxl')
    return df

def plotParamDistributions(df):
    '''
    Purpose: Plot parameter distributions across initial guesses following parameter estimation with 
        training data (Module 2)
    
    Inputs: 
        df: dataframe containing the results of the PEM
        
    Outputs: None
    
    Figures: 
        './FITS Rsq ABOVE 0.99.svg' (plot of training data and simulated data for 
             parameter sets with Rsq > = 0.99)
        'OPTIMIZED PARAMETER DISTRIBUTIONS.svg' (plot of parameter distributions for 
             parameter sets with Rsq > = 0.99)
   '''

    #Only keep rows for which chi2 >
    # df = df[df["Rsq"] >= 0.95]
    all_sims = list(df['Simulation results'])

########################################
# use below to loop through param values in sheet, also loop through each row (for index, row in df.iterrows?)
# then solveall with param values, then plot-- line 543 of Run.py
########################################

    # params_best = []
    # for i in range(0, len(p_all)):
    #     col_name = real_param_labels_all[i] + '*'
    #     val = df[col_name].iloc[0]
    #     params_best.append(val)
    # # print('Params = ', params)
    
    # =============================================================================
    # 1. time series for parameter sets in df
    # ============================================================================
    # fig = plt.figure(figsize = (6.6,2.6))
    # fig.subplots_adjust(wspace=0.1)
    # ax1 = plt.subplot(131)   
    # ax2 = plt.subplot(132)
    # ax3 = plt.subplot(133)
    marker_ = 'o'
    linestyle_ = 'dotted'
    
    exp_1a = exp_data[:6]
    exp_4b = exp_data[6:13]
    exp_4c = exp_data[13:]
    
    err_1a = error[:6]
    err_4b = error[6:13]
    err_4c = error[13:]
    
    t_1a = [0, 24, 48, 72, 96]
    t_4b = [0, 24, 48, 72, 96, 120]
    t_4c = t_4b    
    
    #Plot experimental/training data
    # colors = ['dimgrey', 'black']
    # ax1.errorbar(t_1a, exp_1a[:-1], color = colors[0], marker = marker_, yerr = err_1a[:-1], 
    #              fillstyle = 'none', linestyle = 'none',capsize = 2, label = '1% O2 Training data')
    # ax1.errorbar(t_1a[0], exp_1a[-1], color = colors[1], marker = marker_, yerr = err_1a[-1], 
    #              fillstyle = 'none', linestyle = 'none',capsize = 2, label = '21% O2 Training data')
    
    # ax2.errorbar(t_4b, exp_4b[:-1], color = colors[0], marker = marker_, yerr = err_4b[:-1], 
    #              fillstyle = 'none', linestyle = 'none',capsize = 2, label = '1% O2 Training data')
    # ax2.errorbar(t_4b[0], exp_4b[-1], color = colors[1], marker = marker_, yerr = err_4b[-1], 
    #              fillstyle = 'none', linestyle = 'none',capsize = 2, label = '21% O2 Training data')
    
    # ax3.errorbar(t_4c, exp_4c[:-1], color = colors[0], marker = marker_, yerr = err_4c[:-1], 
    #              fillstyle = 'none', linestyle = 'none',capsize = 2, label = '1% O2 Training data')
    # ax3.errorbar(t_4c[0], exp_4c[-1], color = colors[1], marker = marker_, yerr = err_4c[-1], 
    #          fillstyle = 'none', linestyle = 'none',capsize = 2, label = '21% O2 Training data')
    
    #Plot simulated data for each parameter set in df
    sns.set_palette("Greys", len(all_sims))
    count = 0
    # for sim in all_sims:
    #     count += 1
    #     sim_1a = sim[:6]
    #     sim_4b = sim[6:13]
    #     sim_4c = sim[13:]
        
    #     ax1.plot(t_1a, sim_1a[:-1], marker = None, label = '1% O2 Model fit ' + str(count), 
    #              linestyle = linestyle_)
    #     ax1.plot(t_1a[0], sim_1a[-1], marker = None, label = '21% O2 Model fit ' + str(count), 
    #              linestyle = linestyle_, color = colors[1])
        
    #     ax2.plot(t_4b, sim_4b[:-1], marker = None, label = '1% O2 Model fit ' + str(count), 
    #              linestyle = linestyle_)
    #     ax2.plot(t_4b[0], sim_4b[-1], marker = None, label = '21% O2 Model fit ' + str(count), 
    #              linestyle = linestyle_, color = colors[1])
        
    #     ax3.plot(t_4c, sim_4c[:-1], marker = None, label = '1% O2 Model fit ' + str(count), 
    #              linestyle = linestyle_)
    #     ax3.plot(t_4c[0], sim_4c[-1], marker = None, label = '21% O2 Model fit ' + str(count), 
    #              linestyle = linestyle_, color = colors[1])
 
    # #Set x and y labels and ylim
    # ax1.set_xlabel('Time Post-Plating (hours)')
    # ax1.set_ylabel('Relative DsRE2 Expression')
    # ax1.set_xticks([0, 20, 40, 60, 80, 100])
    # ax1.set_ylim(ax2.get_ylim())
    # ax1.set_title('Simple HBS')
    # ax1.legend()
    # ax1.set_box_aspect(1)
    
    # ax2.set_xlabel('Time Post-Plating (hours)')
    # ax2.set_ylabel('Relative DsRE2 Expression')
    # ax2.set_xticks([0, 20, 40, 60, 80, 100, 120])
    # ax2.set_title('HIF1a Feedback HBS')
    # # ax2.legend()
    # ax2.set_box_aspect(1)
    
    # ax3.set_xlabel('Time Post-Plating (hours)')
    # ax3.set_ylabel('Relative DsRE2 Expression')
    # ax3.set_xticks([0, 20, 40, 60, 80, 100, 120])
    # ax3.set_ylim(ax2.get_ylim())
    # ax3.set_title('HIF2a Feedback HBS')
    # ax3.legend()
    # ax3.set_box_aspect(1)

    # plt.show()
    # plt.savefig('./FITS_Rsq_ABOVE_0.96.svg', bbox_inches="tight")
    
    # =============================================================================
    # 2. parameter distributions for parameter sets in df
    # =============================================================================
    param_labels = []
    for param in real_param_labels_all:
        new_param = param + '*'
        param_labels.append(new_param)

    for label in param_labels:
        new_list = [log10(i) for i in list(df[label])]
        df[label] = new_list
    
    plt.subplots(1,1, figsize=(5.5,4), sharex = True)
    df = pd.melt(df, id_vars=['Rsq'], value_vars=param_labels)
    ax = sns.boxplot(x='variable', y='value', data=df, color = 'gray')
    ax = sns.swarmplot(x='variable', y='value', data=df, color="black")
    ax.set(xlabel='Parameter', ylabel='log(value)')
    
    # plt.show()
    plt.savefig('OPTIMIZED_PARAMETER_DISTRIBUTIONS.svg', dpi = 600)


df = get_df('Opt_Results_chi2_10_pct.xlsx')
# print(df['t_HAF'])
plotParamDistributions(df)