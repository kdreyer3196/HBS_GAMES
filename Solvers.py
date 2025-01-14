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
# MODEL A (base-case)
# =============================================================================

def HBS_1aA(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = HIF2a* protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHH*y6*y1,
            k_txnH*(y4 + y7) - k_dR*y2,
            k_txn - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2*y4 - k_dH1P*y4*y1,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2*y6 - k_bHH*y6*y1,
            k_bHH*y6*y1 - k_dP*y7,
            k_txnBH*(y4 + y7) - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

def HBS_4bA(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = HIF2a* protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHH*y6*y1,
            k_txnH*(y4 + y7) - k_dR*y2,
            k_txn + k_txnBH*(y4 + y7) - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2*y4 - k_dH1P*y4*y1,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2*y6 - k_bHH*y6*y1,
            k_bHH*y6*y1 - k_dP*y7,
            k_txnBH*(y4 + y7) - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

def HBS_4cA(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = HIF2a* protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_bHH*y6*y1,
            k_txnH*(y4 + y7) - k_dR*y2,
            k_txn - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2*y4 - k_dH1P*y4*y1,
            k_txn + k_txnBH*(y4 + y7) - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2*y6 - k_bHH*y6*y1,
            k_bHH*y6*y1 - k_dP*y7,
            k_txnBH*(y4 + y7) - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

# =============================================================================
# MODEL B (add HAF deg in hypoxia to model A)
# =============================================================================

def HBS_1aB(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    if O2 == 138:
        k_dHAF = 0

    else:
        k_dHAF = np.piecewise(t, [t < t_HAF, t >= t_HAF], [k_dHAF, 0])

    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = HIF2a* protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_dHAF*y1 - k_bHH*y6*y1,
            k_txnH*(y4 + y7) - k_dR*y2,
            k_txn - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2*y4 - k_dH1P*y4*y1,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2*y6 - k_bHH*y6*y1,
            k_bHH*y6*y1 - k_dP*y7,
            k_txnBH*(y4 + y7) - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

def HBS_4bB(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    if O2 == 138:
        k_dHAF = 0

    else:
        k_dHAF = np.piecewise(t, [t < t_HAF, t >= t_HAF], [k_dHAF, 0])

    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = HIF2a* protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_dHAF*y1 - k_bHH*y6*y1,
            k_txnH*(y4 + y7) - k_dR*y2,
            k_txn + k_txnBH*(y4 + y7) - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2*y4 - k_dH1P*y4*y1,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2*y6 - k_bHH*y6*y1,
            k_bHH*y6*y1 - k_dP*y7,
            k_txnBH*(y4 + y7) - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

def HBS_4cB(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    if O2 == 138:
        k_dHAF = 0

    else:
        k_dHAF = np.piecewise(t, [t < t_HAF, t >= t_HAF], [k_dHAF, 0])

    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = HIF2a* protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_dHAF*y1 - k_bHH*y6*y1,
            k_txnH*(y4 + y7) - k_dR*y2,
            k_txn - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2*y4 - k_dH1P*y4*y1,
            k_txn + k_txnBH*(y4 + y7) - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2*y6 - k_bHH*y6*y1,
            k_bHH*y6*y1 - k_dP*y7,
            k_txnBH*(y4 + y7) - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt


# =============================================================================
# MODEL C (add SUMOylation to model A)
# =============================================================================

def HBS_1aC(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr
   
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2)*y1*y3,
            k_txn - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2)*y1*y3,
            (k_bHS/O2)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*(y7 + y10) - k_dR*y5,
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4bC(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2)*y1*y3,
            k_txn - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2)*y1*y3,
            (k_bHS/O2)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*(y7 + y10) - k_dR*y5,
            k_txn + k_txnBH*(y7 + y10) - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4cC(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - (k_bHS/O2)*y1*y3,
            k_txn - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2)*y1*y3,
            (k_bHS/O2)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*(y7 + y10) - k_dR*y5,
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2*y7 - k_dH1P*y7*(y1 + y4),
            k_txn + k_txnBH*(y7 + y10) - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt

# =============================================================================
# MODEL D (add HAF deg in hypoxia AND SUMOylation to model A)
# =============================================================================

def HBS_1aD(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    if O2 == 138:
        k_dHAF = 0

    else:
        k_dHAF = np.piecewise(t, [t < t_HAF, t >= t_HAF], [k_dHAF, 0])
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_dHAF*y1 - (k_bHS/O2)*y1*y3,
            k_txn - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2)*y1*y3,
            (k_bHS/O2)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*(y7 + y10) - k_dR*y5,
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4bD(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    if O2 == 138:
        k_dHAF = 0

    else:
        k_dHAF = np.piecewise(t, [t < t_HAF, t >= t_HAF], [k_dHAF, 0])
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_dHAF*y1 - (k_bHS/O2)*y1*y3,
            k_txn - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2)*y1*y3,
            (k_bHS/O2)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*(y7 + y10) - k_dR*y5,
            k_txn + k_txnBH*(y7 + y10) - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2*y7 - k_dH1P*y7*(y1 + y4),
            k_txn - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
    return dydt

def HBS_4cD(y, t, v):

    [[t_HAF, k_txn2, k_dHAF, k_bHS, k_bHH, k_txnH, k_dH1R, k_dH1P, k_dHP, k_txnBH], O2] = v

    #parameters that will be held constant:
    k_txn = 1.0 #U
    k_dR = 2.7 #1/h
    k_tln = 1 #U
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    if O2 == 138:
        k_dHAF = 0

    else:
        k_dHAF = np.piecewise(t, [t < t_HAF, t >= t_HAF], [k_dHAF, 0])
        
    # y holds these state variables: y0 = HAF mRNA, y1 = HAF protein,
    # y2 = SUMO mRNA, y3 = SUMO protein, y4 = SUMO HAF, 
    # y5 = antisense HIF1a RNA, y6 = HIF1a mRNA, y7 = HIF1a protein,
    # y8 = HIF2a mRNA, y9 = HIF2a protein, y10 = HIF2a* protein,
    # y11 = DsRED2 mRNA, y12 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = y

    dydt = [k_txn2 - k_dR*y0,
            k_tln*y0 - k_dP*y1 - k_dHAF*y1 - (k_bHS/O2)*y1*y3,
            k_txn - k_dR*y2,
            k_tln*y2 - k_dP*y3 - (k_bHS/O2)*y1*y3,
            (k_bHS/O2)*y1*y3 - k_dP*y4 - k_bHH*y9*y4,
            k_txnH*(y7 + y10) - k_dR*y5,
            k_txn - k_dR*y6 - k_dH1R*y5*y6,
            k_tln*y6 - k_dP*y7 - k_dHP*O2*y7 - k_dH1P*y7*(y1 + y4),
            k_txn + k_txnBH*(y7 + y10) - k_dR*y8,
            k_tln*y8 - k_dP*y9 - k_dHP*O2*y9 - k_bHH*y9*y4,
            k_bHH*y9*y4 - k_dP*y10,
            k_txnBH*(y7 + y10) - k_dR*y11,
            k_tln*y11 - k_dRep*y12]
        
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
        colors = ['black', 'gray']
            
        if num_states == 9:
            fig, axs = plt.subplots(nrows=3, ncols=3, sharex=False, sharey=False, figsize = (6.5, 6.75))
            fig.subplots_adjust(hspace=0.5)
            fig.subplots_adjust(wspace=0)
            
        elif num_states == 10 or num_states == 11 or num_states == 12:
            fig, axs = plt.subplots(nrows=4, ncols=3, sharex=False, sharey=False, figsize = (6.5, 9))
            fig.subplots_adjust(hspace=0.5)
            fig.subplots_adjust(wspace=0)
            
        elif num_states == 13 or num_states == 14 or num_states == 15:
            fig, axs = plt.subplots(nrows=5, ncols=3, sharex=False, sharey=False, figsize = (6.5, 11.25))
            fig.subplots_adjust(hspace=0.5)
            fig.subplots_adjust(wspace=0)
            
        axs = axs.ravel()
        for i in range(0, num_states):
            axs[i].plot(t_hox, SS_hox[6.6][state_names[i]], color = colors[0], 
                        linestyle = linestyle)
            axs[i].plot(t_hox[0], SS_hox[138][state_names[i]][0], color = colors[1], 
                        linestyle = linestyle)
            axs[i].set_xlabel('Time Post-Plating (hours)')
            axs[i].set_ylabel('1% O2 sim value (a.u.)')
            axs[i].set_title(state_names[i])
            axs[i].set_xticks([0, 20, 40, 60, 80, 100, 120])
            axs[i].ticklabel_format(style='sci', axis='y', scilimits=(-2, 3))
            axs[i].set_box_aspect(1)
            if '1a' in name:
                axs[i].set_xlim([0, 100])

            axs[i].set_ylim(bottom = 0)
            
            max1 = max(SS_hox[6.6][state_names[i]])
            
            # if max1 < 0.5 and max1 >= 0.1:
            #     axs[i].set_ylim(top = max1 + .01 * max1 )
            
            # elif max1 < 0.1:
            #     axs[i].set_ylim(top = max1 + .001 * max1 )
                
            # else:
            axs[i].set_ylim(top = max1 + .1 * max1 )
            
        if num_states == 7:
            axs[-1].axis('off')
            
        elif num_states == 10 or num_states == 13:
            axs[-1].axis('off')
            axs[-2].axis('off')
        
        elif num_states == 11 or num_states == 14:
            axs[-1].axis('off')
                
        fig.suptitle(name + ' States')
        fig_name = name + 'States' + '.svg'
        plt.savefig(fig_name)
        
        print('Timecourses saved as', fig_name)
        
        return t_hox, SS_hox
              
    else:
        return t_hox, SS_hox