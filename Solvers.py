#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 08:39:09 2020

@author: kate
"""

#Package imports
import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# =============================================================================
# CODE TO CALCULATE COST FUNCTIONS
# ============================================================================= 
def calcRsq(data_x, data_y):
    ''' 
    Purpose: Calculate correlation coefficient, Rsq, between 2 datasets
        
    Inputs: 2 lists of floats (of the same length), dataX and dataY 
           
    Output: Rsq value (float) 
    
    '''
    
    #Restructure the data
    x = np.array(data_x)
    y = np.array(data_y)
    x = x.reshape((-1, 1))
    
    #Perform linear regression
    model = LinearRegression()
    model.fit(x,y)
    
    #Calculate Rsq
    Rsq = model.score(x, y)
   
    return Rsq

def calcChi2(exp, sim, std):
    ''' 
    Purpose: 
        Calculate chi2 between 2 datasets with measurement error described by std
        
    Inputs: 
        exp: experimental data (list of floats, length = # datapoints)
        sim: simulated data (list of floats, length = # datapoints)
        std: meaasurement error for exp data (list of floats, length = # datapoints)
           
    Output: 
        chi2: chi2 value (float) 
    
    '''

    #Initialize chi2
    chi2 = 0
    
    #Calculate chi2
    for i, sim_val in enumerate(sim): #for each datapoint
        # err = ((exp[i] - sim_val) / (std[i])) ** 2
        
        err = (exp[i] - sim_val) ** 2
        chi2 = chi2 + err
        
    return chi2

# =============================================================================
# CODE TO DEFINE ODES
# ============================================================================= 

def HBS_1a(y, t, params):

    [k_out, k_prod, k_pMARS, k_dMARS, k_dreg1, k_dreg2, k_dreg3, k_dreg4, O2] = params

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_in = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from hypoxic incubator experiment
    else:
        O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        
    #y holds these state variables: y0= MARS0, y1 = MARS,  y2 = HIF1a mRNA, y3 = HIF1a protein, y4 = HIF2a mRNA,
    #y5 = HIF2a protein, y6 = DsRED2 mRNA, y7 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7 = y

    dydt = [k_in + k_prod*(y3 + y5) - k_pMARS*y0 - k_out*y0,
            k_pMARS*y0 +  - k_dMARS*y1,
            k_txn - k_dR*y2 - k_dreg1*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dreg2*y1*O2_rate*y3,
            k_txn - k_dR*y4 + k_dreg3*y1,
            k_tln*y4 - k_dP*y5 - k_dreg4*y1*O2_rate*y5,
            
            k_txnBh*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_drep*y7]

    return dydt

def HBS_4b(y, t, params):

    [k_out, k_prod, k_pMARS, k_dMARS, k_dreg1, k_dreg2, k_tln2, k_dreg3, k_dreg4, O2] = params

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr
    
    k_in = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    
        #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
    else:
        #Exponential decrease in pO2 from hypoxic incubator experiment
        O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        
    #y holds these state variables: y0 = MARS0, y1= MARS,  y2 = HIF1a mRNA, y3 = synthetic HIF1a mRNA, y4 = HIF1a protein,
    #y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_in + k_prod*(y4 + y6) - k_pMARS*y0 - k_out*y0,
            k_pMARS*y0 +  - k_dMARS*y1,
            k_txn - k_dR*y2 - k_dreg1*y1*y2,
            k_txnBh*(y4 + y6) - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dreg2*y1*O2_rate*y4,
            k_txn - k_dR*y5 + k_dreg3*y1,
            k_tln*y5 - k_dP*y6 - k_dreg4*y1*O2_rate*y6,
            
            k_txnBh*(y4 + y6) - k_dR*y7,
            k_tln*y7 - k_drep*y8]

    return dydt

def HBS_4c(y, t, params):

    [k_out, k_prod, k_pMARS, k_dMARS, k_dreg1, k_dreg2, k_tln2, k_dreg3, k_dreg4, O2] = params

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr
    
    k_in = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

        #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
    else:
        #Exponential decrease in pO2 from hypoxic incubator experiment
        O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        
    #y holds these state variables: Y0 = MARS0, y1 = MARS,  y2 = HIF1a mRNA, y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = synthetic HIF2a mRNA,
    #y6 = HIF2a protein, y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_in + k_prod*(y3 + y6) - k_pMARS*y0 - k_out*y0,
            k_pMARS*y0 +  - k_dMARS*y1,
            k_txn - k_dR*y2 - k_dreg1*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dreg2*y1*O2_rate*y3,
            k_txn - k_dR*y4 + k_dreg3*y1,
            k_txnBh*(y3 + y6) - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dreg4*y1*O2_rate*y6,

            k_txnBh*(y3 + y6) - k_dR*y7,
            k_tln*y7 - k_drep*y8]

    return dydt


# =============================================================================
# CODE TO SOLVE ONE HBS in Normoxia and Hypoxia (START HERE)
# ============================================================================= 
def solveSingle(args): 
    
    ''' 
    Purpose: Given a set of arguments (parameters, initial conditions), 
    solve the model for a single datapoint (for a single set of component doses)
        
    Inputs: 
        args (a list of arguments defining the condition) defined as
        args = [index, v, tp1_, tp2_,  output, model]
           
    Output: 
        if output == 'save internal states'
            t_1 (list of timepoints before ligand addition)
            t_2 (list of timepoints after ligand addition)
            solution_before (2D array, dynamic trajectory of all states before ligand addition)
            solution_on (2D array, dynamic trajectory of all states after ligand addition)
            state_labels (list of labels for each state variable in the model)
        
        OR 
        
        if output != 'save internal states'
            solution_on[-1, -1] (float, the value of the reporter state at the final timepoint)
    
    Figures: 
        if output == 'timecourse':
            './TIMECOURSE.svg' (plots of dynamic trajectories of each state variable)
    
    '''
    
    [index, v, tp1_, tp2_, output, model] = args

    #We only consider a single model structure in this work. 
    #To change the model structure (without changing the free params or training data)
    #Add an elif statement to this section of the code with the appropirate information. 
    
    num_states = 8
    output_state = 7
    state_labels = ['A mRNA', 'A protein', 'B mRNA', 'B protein', 
                    'Ligand', 'TF', 'Reporter mRNA', 'Reporter protein']
    
    if model in ['model A', 'model B']:
        ode_model = model_AB
        
    elif model in ['model C', 'model D']:
        ode_model = model_CD

    else: 
        print('Error: Model name does not exist.')

    #Define ligand dose
    dose_ligand = v[0][2]

    #Set time 
    numpoints = 100 #number of timepoints
    t_1 = [tp1_ * float(i) / (numpoints - 1) for i in range(numpoints)]  #range of timepoints
    
    #Set initial conditions to 0
    y0_before = [0] * num_states
    
    #1. Solve the ODEs for the time up until ligand addition
    solution_before = odeint(ode_model, y0_before, t_1, args=(v,), mxstep=50000000)
    
    #2. Solve equations at and after Ligand addition timepoint
    #set ligand initial condition to dose_ligand
    solution_before[-1, 4] =  dose_ligand
    
    #set initial conditions to the final timepoint values of 
    #each state variable from the "before" condition
    y0_on = solution_before[-1, :]

    #Set time
    numpoints = 25
    t_2_ = [tp2_ * float(i) / (numpoints - 1) for i in range(numpoints)]  # hours
  
    #Solve equations after Ligand addition 
    solution_on = odeint(ode_model, y0_on, t_2_, args=(v,), mxstep=500000000)
    
    #Define t_2 for plotting
    t_2 = [i + tp1_ for i in t_2_]  

    if output == 'timecourse':
        
        #Plot timecourse
        fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=False, figsize = (8, 4))
        fig.subplots_adjust(hspace=.5)
        fig.subplots_adjust(wspace=0.3)
        
        axs = axs.ravel()
        for i in range(0, num_states):
            axs[i].plot(t_1, solution_before[:,i], color = 'black', linestyle = 'dashed')
            axs[i].plot(t_2, solution_on[:,i], color = 'black')
            axs[i].set_xlabel('Time (hours)', fontsize = 8)
            axs[i].set_ylabel('Simulation value (a.u.)', fontsize = 8)
            axs[i].set_title(state_labels[i], fontweight = 'bold', fontsize = 10)
            
            max1 = max(solution_before[:,i])
            max2 = max(solution_on[:,i])
            axs[i].set_ylim(top = max(max1, max2) + .1 * max(max1, max2) )
        plt.savefig('./TIMECOURSES.svg')
        
        return 'Timecourses saved as TIMECOURSES.svg'
              
    elif output == 'save internal states':
        return t_1, t_2, solution_before, solution_on, state_labels
    else:
        return solution_on[-1, output_state] 


    
