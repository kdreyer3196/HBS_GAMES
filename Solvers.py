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
# MODEL 4C
# =============================================================================

def HBS_1a4C(y, t, v):

    [[k_bHS, k_rbHS, k_txnH, k_txnHAF, k_dH1R, k_dH1P, k_dHP, k_txnBH, k_aH2P, k_iH2P, k_txnb1, k_txnb2], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_txn = 1 #au/h
    # k_txnH = k_txn

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, #y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein,  y10 = DsRED2 mRNA
    # y11 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3,
            (k_bHS/O2_rate)*y1*y3 - k_dP*y4,
            k_txnH*y7 + k_txnH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1 - k_dR*y5,
            
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9,
            k_txnBH*y7 + k_txnBH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1- k_dR*y10,
            
            k_tln*y10 - k_dRep*y11]
        
    return dydt

def HBS_4b4C(y, t, v):

    [[k_bHS, k_rbHS, k_txnH, k_txnHAF, k_dH1R, k_dH1P, k_dHP, k_txnBH, k_aH2P, k_iH2P, k_txnb1, k_txnb2], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_txn = 1 #au/h
    # k_txnH = k_txn

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, #y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein,  y10 = DsRED2 mRNA
    # y11 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3,
            (k_bHS/O2_rate)*y1*y3 - k_dP*y4,
            k_txnH*y7 + k_txnH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1 - k_dR*y5,
            
            k_txn + k_txnBH*y7 + k_txnBH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1 - k_dR*y6 - k_dH1R*y5*y6,
            
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9,
            k_txnBH*y7 + k_txnBH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1- k_dR*y10,
            
            k_tln*y10 - k_dRep*y11]
        
    return dydt

def HBS_4c4C(y, t, v):

    [[k_bHS, k_rbHS, k_txnH, k_txnHAF, k_dH1R, k_dH1P, k_dHP, k_txnBH, k_aH2P, k_iH2P, k_txnb1, k_txnb2], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_txn = 1 #au/h
    # k_txnH = k_txn

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, #y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein,  y10 = DsRED2 mRNA
    # y11 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3,
            (k_bHS/O2_rate)*y1*y3 - k_dP*y4,
            k_txnH*y7 + k_txnH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1 - k_dR*y5,
            
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txn + k_txnBH*y7 + k_txnBH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1 - k_dR*y8,
            
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9,
            k_txnBH*y7 + k_txnBH*(k_aH2P*y9*y4/(1 + k_aH2P*y9*y4 + k_iH2P*y9))
            + k_txnHAF*y1- k_dR*y10,
            
            k_tln*y10 - k_dRep*y11]
        
    return dydt

# =============================================================================
# MODEL 4D
# =============================================================================

def HBS_1a4D(y, t, v):

    [[k_bHS, k_rbHS, k_txnH, k_dH1R, k_dH1P, k_dHP, k_bHH, k_txnBH, k_aH2P, k_iH2P, k_txnb1, k_txnb2, k_txn], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    # k_txn = 1 #au/h
    # k_txnH = k_txn

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3,
            (k_bHS/O2_rate)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*y7 + k_txnH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y5,
            
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            
            k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y11,
            
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4b4D(y, t, v):

    [[k_bHS, k_rbHS, k_txnH, k_dH1R, k_dH1P, k_dHP, k_bHH, k_txnBH, k_aH2P, k_iH2P, k_txnb1, k_txnb2, k_txn], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    # k_txn = 1 #au/h
    # k_txnH = k_txn

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3,
            (k_bHS/O2_rate)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*y7 + k_txnH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y5,
            
            k_txn + k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y6 - k_dH1R*y5*y6,
            
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            
            
            k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y11,
            
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4c4D(y, t, v):
    
    [[k_bHS, k_rbHS, k_txnH, k_dH1R, k_dH1P, k_dHP, k_bHH, k_txnBH, k_aH2P, k_iH2P, k_txnb1, k_txnb2, k_txn], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    # k_txn = 1 #au/h
    # k_txnH = k_txn

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3,
            (k_bHS/O2_rate)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*y7 + k_txnH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y5,
            
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txn + k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y8,
            
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            
            k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y11,
            
            k_tln*y11 - k_dRep*y12]
        
    return dydt

# =============================================================================
# MODEL 4E
# =============================================================================

def HBS_1a4E(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            (k_bHS/O2_rate)*y1*y3 - k_rbHS*y4 - k_dP*y4 - k_bHH*y9*y4
            + k_rbHH*y10,
            
            k_txnH*y7 + k_txnH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y5,
            
            k_txnb - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txnb3 - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4 + k_rbHH*y10,
            k_bHH*y9*y4 - k_rbHH*y10 - k_dP*y10,
            k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y11,
            
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4b4E(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            (k_bHS/O2_rate)*y1*y3 - k_rbHS*y4 - k_dP*y4 - k_bHH*y9*y4
            + k_rbHH*y10,
            
            k_txnH*y7 + k_txnH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y5,
            
            k_txnb + k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y6 - k_dH1R*y5*y6,
            
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txnb3 - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4 + k_rbHH*y10,
            k_bHH*y9*y4 - k_rbHH*y10 - k_dP*y10,
            k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y11,
            
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4c4E(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            (k_bHS/O2_rate)*y1*y3 - k_rbHS*y4 - k_dP*y4 - k_bHH*y9*y4
            + k_rbHH*y10,
            
            k_txnH*y7 + k_txnH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y5,
            
            k_txnb - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txnb3 + k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y8,
            
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4 + k_rbHH*y10,
            k_bHH*y9*y4 - k_rbHH*y10 - k_dP*y10,
            k_txnBH*y7 + k_txnBH*(k_aH2P*y10/(1 + k_aH2P*y10 + k_iH2P*y9))
            - k_dR*y11,
            
            k_tln*y11 - k_dRep*y12]
        
    return dydt

# =============================================================================
# MODEL 4F
# =============================================================================

def HBS_1a4F(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            (k_bHS/O2_rate)*y1*y3 - k_rbHS*y4 - k_dP*y4 - k_bHH*y9*y4
            + k_rbHH*y10,
            
            k_txnH*(y7 + y10) - k_dR*y5 - k_dH1R*y5*y6,
            k_txnb - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txnb3 - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4 + k_rbHH*y10,
            k_bHH*y9*y4 - k_rbHH*y10 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt


def HBS_4b4F(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            (k_bHS/O2_rate)*y1*y3 - k_rbHS*y4 - k_dP*y4 - k_bHH*y9*y4
            + k_rbHH*y10,
            
            k_txnH*(y7 + y10) - k_dR*y5 - k_dH1R*y5*y6,
            k_txnb + k_txnBH*(y7 + y10) - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txnb3 - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4 + k_rbHH*y10,
            k_bHH*y9*y4 - k_rbHH*y10 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4c4F(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2_rate)*y1*y3 + k_rbHS*y4,
            (k_bHS/O2_rate)*y1*y3 - k_rbHS*y4 - k_dP*y4 - k_bHH*y9*y4
            + k_rbHH*y10,
            
            k_txnH*(y7 + y10) - k_dR*y5 - k_dH1R*y5*y6,
            k_txnb - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2_rate*y7 - k_dH1P*y7*(y1 + y4),
            k_txnb3 + k_txnBH*(y7 + y10) - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2_rate*y9 - k_bHH*y9*y4 + k_rbHH*y10,
            k_bHH*y9*y4 - k_rbHH*y10 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt

# =============================================================================
# MODEL 5
# =============================================================================

def HBS_1a5(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_m, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO* protein,
    # y5 = SUMO HAF, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHS*y1*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_m/O2_rate)*y3,
            (k_m/O2_rate)*y3 - k_dP*y4 - k_bHS*y1*y4 + k_rbHS*y5,
            k_bHS*y1*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6 - k_dH1R*y6*y7,
            k_txnb - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y1 + y5),
            k_txnb3 - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

def HBS_4b5(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_m, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO* protein,
    # y5 = SUMO HAF, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHS*y1*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_m/O2_rate)*y3,
            (k_m/O2_rate)*y3 - k_dP*y4 - k_bHS*y1*y4 + k_rbHS*y5,
            k_bHS*y1*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6 - k_dH1R*y6*y7,
            k_txnb + k_txnBH*(y8 + y11) - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y1 + y5),
            k_txnb3 - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

def HBS_4c5(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_m, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO* protein,
    # y5 = SUMO HAF, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHS*y1*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_m/O2_rate)*y3,
            (k_m/O2_rate)*y3 - k_dP*y4 - k_bHS*y1*y4 + k_rbHS*y5,
            k_bHS*y1*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6 - k_dH1R*y6*y7,
            k_txnb - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y1 + y5),
            k_txnb3 + k_txnBH*(y8 + y11) - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

# =============================================================================
# MODEL 5A
# =============================================================================

def HBS_1a5A(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_m, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO* protein,
    # y5 = HAF_s, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHS*y1*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_m/O2_rate)*y3,
            (k_m/O2_rate)*y3 - k_dP*y4 - k_bHS*y1*y4 + k_rbHS*y5,
            k_bHS*y1*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6,
            k_txnb - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y1 + y5),
            k_txnb3 - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

def HBS_4b5A(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_m, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO* protein,
    # y5 = SUMO HAF, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHS*y1*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_m/O2_rate)*y3,
            (k_m/O2_rate)*y3 - k_dP*y4 - k_bHS*y1*y4 + k_rbHS*y5,
            k_bHS*y1*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6,
            k_txnb + k_txnBH*(y8 + y11) - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y1 + y5),
            k_txnb3 - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

def HBS_4c5A(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_m, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO* protein,
    # y5 = SUMO HAF, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHS*y1*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_m/O2_rate)*y3,
            (k_m/O2_rate)*y3 - k_dP*y4 - k_bHS*y1*y4 + k_rbHS*y5,
            k_bHS*y1*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6,
            k_txnb - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y1 + y5),
            k_txnb3 + k_txnBH*(y8 + y11) - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

# =============================================================================
# MODEL 6
# =============================================================================

def HBS_1a6(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_p, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
        k_p = 0

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
        k_p = np.piecewise(t, [t < 44, t >= 44], [0, k_p])
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = HAF* protein, y3 = SUMO mRNA, y4 = SUMO protein,
    # y5 = HAF_s, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_p*y1,
            k_p*y1 - k_dP*y2 - k_bHS*y2*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y3,
            k_tln*y3 - k_dP*y4,
            k_bHS*y2*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6,
            k_txnb - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y2 + y5),
            k_txnb3 - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

def HBS_4b6(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_p, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
        k_p = 0

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
        k_p = np.piecewise(t, [t < 44, t >= 44], [0, k_p])
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = HAF* protein, y3 = SUMO mRNA, y4 = SUMO protein,
    # y5 = HAF_s, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_p*y1,
            k_p*y1 - k_dP*y2 - k_bHS*y2*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y3,
            k_tln*y3 - k_dP*y4,
            k_bHS*y2*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6,
            k_txnb + k_txnBH*(y8 + y11) - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y2 + y5),
            k_txnb3 - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
    return dydt

def HBS_4c6(y, t, v):

    [[k_txnb1, k_bHS, k_rbHS, k_txnb2, k_p, k_bHH, k_rbHH, k_txnH, k_txnb3, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    # k_txnBH = 1.0
    k_txnb = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
        k_p = 0

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
        k_p = np.piecewise(t, [t < 44, t >= 44], [0, k_p])
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = HAF* protein, y3 = SUMO mRNA, y4 = SUMO protein,
    # y5 = HAF_s, y6 = antisense HIF1a RNA, y7 = HIF1a mRNA,
    # y8 = HIF1a protein, y9 = HIF2a mRNA, y10 = HIF2a protein,
    # y11 = HIF2a* protein, y12 = DsRED2 mRNA, y13 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13 = y

    dydt = [k_txnb1 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_p*y1,
            k_p*y1 - k_dP*y2 - k_bHS*y2*y4 + k_rbHS*y5,
            k_txnb2 - k_dR*y3,
            k_tln*y3 - k_dP*y4,
            k_bHS*y2*y4 - k_rbHS*y5 - k_dP*y5 - k_bHH*y10*y5
            + k_rbHH*y11,
            
            k_txnH*(y8 + y11) - k_dR*y6,
            k_txnb - k_dR*y7 - k_dH1R*y6*y7,
            k_tln*y7 - k_dP*y8 - k_dHP*O2_rate*y8 - k_dH1P*y8*(y2 + y5),
            k_txnb3 + k_txnBH*(y8 + y11) - k_dR*y9,
            k_tln*y9 - k_dP*y10 - k_dHP*O2_rate*y10 - k_bHH*y10*y5 + k_rbHH*y11,
            k_bHH*y10*y5 - k_rbHH*y11 - k_dP*y11,
            k_txnBH*(y8 + y11) - k_dR*y12,
            k_tln*y12 - k_dRep*y13]
        
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
        [ODEs, v, name, num_states, state_names, output, O2_range, t_hox] = args           
 
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
        if output == 'all states':
            './All_States_model_X.svg' (plots of dynamic trajectories of each state variable)
    
    '''
    
    [ODEs, v, name, num_states, state_names, output, O2_range, t_hox] = args
 
    y0_nox = np.zeros(num_states)
    t_ss = np.arange(0, 500, 1) #hours
    v[1] = 138
    
    solution_nox = spi.odeint(ODEs, y0_nox, t_ss, args=(v,))
    
    ###create dictionary to store SS normoxia solutions
    SS_nox = {}
    
    for i in range(0,num_states):
        name_st = state_names[i]
        
        SS_nox[name_st] = solution_nox[:,i]

    ###hypoxia simulation###
    y0_hox = []
    
    for name_st in state_names:
        y0_hox.append(SS_nox[name_st][-1])
            
    ###create dictionary to store SS normoxia solutions
    SS_hox = {}
    for O2_val in O2_range:
        v[1] = O2_val
        SS_hox[O2_val] = {}
        solution_hox = spi.odeint(ODEs, y0_hox, t_hox, args=(v,))
        
        for i in range(0,num_states):
            name_st = state_names[i]
            
            SS_hox[O2_val][name_st] = solution_hox[:,i]
#############################################################################################

    if output == 'all states':
        
        #Plot timecourse
        linestyle = 'dotted'
        sky_blue = [i/255 for i in [86, 180, 233]]
        colors = [sky_blue, 'gray']
        
        if num_states == 7:
            fig, axs = plt.subplots(nrows=2, ncols=4, sharex=False, sharey=False, figsize = (8, 4))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
        
        elif num_states == 8:
            fig, axs = plt.subplots(nrows=2, ncols=4, sharex=False, sharey=False, figsize = (8, 4))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
            
        elif num_states == 9:
            fig, axs = plt.subplots(nrows=3, ncols=3, sharex=False, sharey=False, figsize = (6, 6))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
            
        elif num_states == 10:
            fig, axs = plt.subplots(nrows=2, ncols=5, sharex=False, sharey=False, figsize = (10, 4))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
            
        elif num_states == 12:
            fig, axs = plt.subplots(nrows=3, ncols=4, sharex=False, sharey=False, figsize = (8, 6))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
            
        elif num_states == 13 or num_states == 14:
            fig, axs = plt.subplots(nrows=3, ncols=5, sharex=False, sharey=False, figsize = (10, 6))
            fig.subplots_adjust(hspace=.5)
            fig.subplots_adjust(wspace=0.3)
            
        axs = axs.ravel()
        for i in range(0, num_states):
            axs[i].plot(t_hox, SS_hox[7.6][state_names[i]], color = colors[0], 
                        linestyle = linestyle)
            axs[i].plot(t_hox[0], SS_hox[138][state_names[i]][0], color = colors[1], 
                        linestyle = linestyle)
            axs[i].set_xlabel('Time Post-Plating (hours)')
            axs[i].set_ylabel('1% O2 sim value (a.u.)')
            axs[i].set_title(state_names[i])
            axs[i].set_ylim(bottom = 0)
            
            max1 = max(SS_hox[7.6][state_names[i]])
            
            if max1 < 0.5 and max1 >= 0.1:
                axs[i].set_ylim(top = max1 + .01 * max1 )
            
            elif max1 < 0.1:
                axs[i].set_ylim(top = max1 + .001 * max1 )
                
            else:
                axs[i].set_ylim(top = max1 + .1 * max1 )
            
        if num_states == 7:
            axs[-1].axis('off')
            
        elif num_states == 13:
            axs[-1].axis('off')
            axs[-2].axis('off')
        
        elif num_states == 14:
            axs[-1].axis('off')
                
        fig.suptitle(name + ' States')
        fig_name = name + 'States' + '.svg'
        plt.savefig(fig_name)
        
        print('Timecourses saved as', fig_name)
        
        return t_hox, SS_hox
              
    else:
        return t_hox, SS_hox