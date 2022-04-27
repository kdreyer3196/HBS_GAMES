#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 08:39:09 2020

@author: kate
"""

#Package imports
import math
import numpy as np
import scipy.integrate as spi
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

# =============================================================================
# MODEL 0
# =============================================================================

def HBS_1a(y, t, v):
 
    [[k_prod, k_pMARS, k_out, k_dMARS, k_dreg1, k_tln2, k_dreg2, k_dreg3, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_in = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dreg4 = k_dreg2*deg_ratio
    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
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

def HBS_4b(y, t, v):

    [[k_prod, k_pMARS, k_out, k_dMARS, k_dreg1, k_tln2, k_dreg2, k_dreg3, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr
    
    k_in = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    
    k_dreg4 = k_dreg2*deg_ratio
    
        #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
    else:
        #Exponential decrease in pO2 from fit to diffusion equation
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
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

def HBS_4c(y, t, v):

    [[k_prod, k_pMARS, k_out, k_dMARS, k_dreg1, k_tln2, k_dreg2, k_dreg3, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr
    
    k_in = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    
    k_dreg4 = k_dreg2*deg_ratio

        #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
    else:
        #Exponential decrease in pO2 from fit to diffusion equation
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
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
# MODEL 1
# =============================================================================

def HBS_1a1(y, t, v):
 
    [[k_pMARS, k_dMARS, k_dreg1, k_dO2, k_dreg2, k_act], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0= MARS, y1 = HIF1a mRNA, y2 = HIF1a
    # protein, y3 = HIF2a mRNA, y4 = HIF2a protein, y5 = DsRED2 mRNA,
    # y6 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6 = y

    dydt = [k_basal + k_pMARS*(y2 + y4) - k_dMARS*y0,
            k_txn - k_dR*y1 - k_dreg1*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dO2*O2_rate*y2 - k_dreg2*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dO2*O2_rate*y4,
            k_txnBh*y2 + k_txnBh*y4*(k_act*y0/(1 + k_act*y0)) - k_dR*y5,
            k_tln*y5 - k_drep*y6]
        
    return dydt

def HBS_4b1(y, t, v):
 
    [[k_pMARS, k_dMARS, k_dreg1, k_dO2, k_dreg2, k_act], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0= MARS, y1 = HIF1a mRNA, y2 = HIF1a
    # protein, y3 = HIF2a mRNA, y4 = HIF2a protein, y5 = DsRED2 mRNA,
    # y6 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6 = y

    dydt = [k_basal + k_pMARS*(y2 + y4) - k_dMARS*y0,
            k_txn + k_txnBh*y2 + k_txnBh*y4*(k_act*y0/(1 + k_act*y0)) - k_dR*y1
            - k_dreg1*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dO2*O2_rate*y2 - k_dreg2*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dO2*O2_rate*y4,
            k_txnBh*y2 + k_txnBh*y4*(k_act*y0/(1 + k_act*y0)) - k_dR*y5,
            k_tln*y5 - k_drep*y6]
        
    return dydt

def HBS_4c1(y, t, v):
 
    [[k_pMARS, k_dMARS, k_dreg1, k_dO2, k_dreg2, k_act], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0= MARS, y1 = HIF1a mRNA, y2 = HIF1a
    # protein, y3 = HIF2a mRNA, y4 = HIF2a protein, y5 = DsRED2 mRNA,
    # y6 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6 = y

    dydt = [k_basal + k_pMARS*(y2 + y4) - k_dMARS*y0,
            k_txn - k_dR*y1 - k_dreg1*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dO2*O2_rate*y2 - k_dreg2*y0*y2,
            k_txn + k_txnBh*y2 + k_txnBh*y4*(k_act*y0/(1 + k_act*y0))
            - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dO2*O2_rate*y4,
            k_txnBh*y2 + k_txnBh*y4*(k_act*y0/(1 + k_act*y0)) - k_dR*y5,
            k_tln*y5 - k_drep*y6]
        
    return dydt

# =============================================================================
# MODEL 2
# =============================================================================

def HBS_1a2(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dHP, k_dH1P, k_aH2P], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = MARS0, y1 = MARS1, y2 = HIF1a mRNA,
    # y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = HIF2a protein, y6 = DsRED2 mRNA,
    # y7 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7 = y

    dydt = [k_basal + k_pM0*(y3 + y5) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBh*y3 + k_txnBh*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_drep*y7]
        
    return dydt

def HBS_4b2(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dHP, k_dH1P, k_aH2P], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = MARS0, y1 = MARS1, y2 = HIF1a mRNA,
    # y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = HIF2a protein, y6 = DsRED2 mRNA,
    # y7 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7 = y

    dydt = [k_basal + k_pM0*(y3 + y5) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 - k_dM1*y1,
            k_txn + k_txnBh*y3 + k_txnBh*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) -
            k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBh*y3 + k_txnBh*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_drep*y7]
        
    return dydt

def HBS_4c2(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dHP, k_dH1P, k_aH2P], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = MARS0, y1 = MARS1, y2 = HIF1a mRNA,
    # y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = HIF2a protein, y6 = DsRED2 mRNA,
    # y7 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7 = y

    dydt = [k_basal + k_pM0*(y3 + y5) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn  + k_txnBh*y3 + k_txnBh*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) -
            k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBh*y3 + k_txnBh*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_drep*y7]
        
    return dydt

# =============================================================================
# MODEL 3
# =============================================================================

def HBS_1a3(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dHP, k_dH1P, k_aH2P, k_tln2], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = MARS0, y1 = MARS1, y2 = HIF1a mRNA,
    # y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = HIF2a protein, y6 = DsRED2 mRNA,
    # y7 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7 = y

    dydt = [k_basal + k_pM0*(y3 + y5) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBh*y3 + k_txnBh*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_drep*y7]
        
    return dydt

def HBS_4b3(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dHP, k_dH1P, k_aH2P, k_tln2], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = MARS0, y1 = MARS1, y2 = HIF1a mRNA,
    # y3 = sHIF1a mRNA, y4 = HIF1a protein, y5 = HIF2a mRNA, y6 = HIF2a protein,
    # y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_basal + k_pM0*(y4 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_txnBh*y4 + k_txnBh*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dHP*O2_rate*y4 - k_dH1P*y1*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2_rate*y6,
            k_txnBh*y4 + k_txnBh*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y7,
            k_tln*y7 - k_drep*y8]
        
    return dydt

def HBS_4c3(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dHP, k_dH1P, k_aH2P, k_tln2], O2] = v

    #parameters that will be held constant:
    k_txnBh = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_drep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = MARS0, y1 = MARS1, y2 = HIF1a mRNA,
    # y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = sHIF2a mRNA, y6 = HIF2a protein,
    # y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_basal + k_pM0*(y3 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_txnBh*y3 + k_txnBh*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dHP*O2_rate*y6,
            
            k_txnBh*y3 + k_txnBh*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y7,
            k_tln*y7 - k_drep*y8]
    
    return dydt


# =============================================================================
# CODE TO SOLVE ONE HBS in Normoxia and Hypoxia 
# ============================================================================= 
def solveSingle(args): 
    
    ''' 
    Purpose: Given a set of arguments (parameters, initial conditions), 
    solve the model for a single datapoint (for a single set of component doses)
        
    Inputs: 
        args (a list of arguments defining the condition) defined as
        args = [model, v, num_states, state_names, O2_range, output]
           
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
    
    [ODEs, v, num_states, state_names, output, O2_range, t_hox] = args

    y0_nox = np.zeros(num_states)
    t_ss = np.arange(0, 500, 1) #hours
    v[1] = 138
    
    solution_nox = spi.odeint(ODEs, y0_nox, t_ss, args=(v,))
    
    ###create dictionary to store SS normoxia solutions
    SS_nox = {}
    
    for i in range(0,num_states):
        name = state_names[i]
        
        SS_nox[name] = solution_nox[:,i]

    ###hypoxia simulation###
    y0_hox = []
    
    for name in state_names:
        y0_hox.append(SS_nox[name][-1])
            
    ###create dictionary to store SS normoxia solutions
    SS_hox = {}
    for O2_val in O2_range:
        v[1] = O2_val
        SS_hox[O2_val] = {}
        solution_hox = spi.odeint(ODEs, y0_hox, t_hox, args=(v,))
        
        for i in range(0,num_states):
            name = state_names[i]
            
            SS_hox[O2_val][name] = solution_hox[:,i]
#############################################################################################

    if output == 'timecourse':
        
        #Plot timecourse
        if num_states == 8:
            fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=False, figsize = (8, 4))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
            
        if num_states == 9:
            fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=False, figsize = (6, 6))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
            
        axs = axs.ravel()
        for i in range(0, num_states):
            axs[i].plot(t_hox, SS_hox[7.6][state_names[i]], color = 'dimgrey', 
                        label = '1% O2')
            axs[i].plot(t_hox, SS_hox[138][state_names[i]], color = 'black', 
                        label = '21% O2')
            axs[i].set_xlabel('Time (hours)', fontsize = 8)
            axs[i].set_ylabel('Simulation value (a.u.)', fontsize = 8)
            axs[i].set_title(state_names[i], fontweight = 'bold', fontsize = 10)
            
            max1 = max(SS_hox[7.6][state_names[i]])
            axs[i].set_ylim(top = max1 + .1 * max1 )
        plt.savefig('./TIMECOURSES.svg')
        
        return 'Timecourses saved as TIMECOURSES.svg'
              
    else:
        return t_hox, SS_hox


    
