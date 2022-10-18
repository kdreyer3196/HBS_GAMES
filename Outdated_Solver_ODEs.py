# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:54:50 2022

@author: Katie_Dreyer
"""
import numpy as np

# =============================================================================
# MODEL 0
# =============================================================================

def HBS_1a0(y, t, v):
     
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio
    
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

    dydt = [k_basal + k_pM0*(y3 + y5) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dH1P*y1*O2_rate*y3,
            k_txn - k_dR*y4 + k_pH2R*y1,
            k_tln*y4 - k_dP*y5 - k_dH2P*y1*O2_rate*y5,
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]

    return dydt

def HBS_4b0(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio
    
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

    dydt = [k_basal + k_pM0*(y4 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_txnBH*(y4 + y6) - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dH1P*y1*O2_rate*y4,
            k_txn - k_dR*y5 + k_pH2R*y1,
            k_tln*y5 - k_dP*y6 - k_dH2P*y1*O2_rate*y6, 
            k_txnBH*(y4 + y6) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]

    return dydt

def HBS_4c0(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio

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

    dydt = [k_basal + k_pM0*(y3 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dH1P*y1*O2_rate*y3,
            k_txn - k_dR*y4 + k_pH2R*y1,
            k_txnBH*(y3 + y6) - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dH2P*y1*O2_rate*y6,
            k_txnBH*(y3 + y6) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]

    return dydt

# =============================================================================
# MODEL 0A
# =============================================================================

def HBS_1a0A(y, t, v):
     
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio
    
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

    dydt = [k_basal + k_pM0*(y3 + y5) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dH1P*y1*O2_rate*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dH2P*y1*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)),
            k_tln*y6 - k_dRep*y7]

    return dydt

def HBS_4b0A(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio
    
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

    dydt = [k_basal + k_pM0*(y4 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_txnBH*y4 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dH1P*y1*O2_rate*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dH2P*y1*O2_rate*y6, 
            k_txnBH*y4 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]

    return dydt

def HBS_4c0A(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio

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

    dydt = [k_basal + k_pM0*(y3 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dH1P*y1*O2_rate*y3,
            k_txn - k_dR*y4,
            k_txnBH*y3 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dH2P*y1*O2_rate*y6,
            k_txnBH*y3 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]

    return dydt

# =============================================================================
# MODEL 0B
# =============================================================================

def HBS_1a0B(y, t, v):
     
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio
    
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

    dydt = [k_basal + k_pM0*(y3 + y5) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dH1P*y1*O2_rate*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dH2P*y1*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*y5*k_aH2P*y1,
            k_tln*y6 - k_dRep*y7]

    return dydt

def HBS_4b0B(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio
    
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

    dydt = [k_basal + k_pM0*(y4 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_txnBH*y4 + k_txnBH*y6*k_aH2P*y1 - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dH1P*y1*O2_rate*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dH2P*y1*O2_rate*y6, 
            k_txnBH*y4 + k_txnBH*y6*k_aH2P*y1 - k_dR*y7,
            k_tln*y7 - k_dRep*y8]

    return dydt

def HBS_4c0B(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h

    k_dH2P = k_dH1P*deg_ratio

        #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2
    else:
        #Exponential decrease in pO2 from fit to diffusion equation
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    #y holds these state variables: y0 = MARS0, y1 = MARS,  y2 = HIF1a mRNA, y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = synthetic HIF2a mRNA,
    #y6 = HIF2a protein, y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_basal + k_pM0*(y3 + y6) - k_pM1*y0 - k_dM0*y0,
            k_pM1*y0 +  - k_dM1*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dH1P*y1*O2_rate*y3,
            k_txn - k_dR*y4,
            k_txnBH*y3 + k_txnBH*y6*k_aH2P*y1 - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dH2P*y1*O2_rate*y6,
            k_txnBH*y3 + k_txnBH*y6*k_aH2P*y1 - k_dR*y7,
            k_tln*y7 - k_dRep*y8]

    return dydt

# =============================================================================
# MODEL 1
# =============================================================================

def HBS_1a1(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*y2 + k_txnBH*y4*(k_aH2P*y0/(1 + k_aH2P*y0)) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4b1(y, t, v):
    
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn + k_txnBH*y2 + k_txnBH*y4*(k_aH2P*y0/(1 + k_aH2P*y0)) - k_dR*y1
            - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*y2 + k_txnBH*y4*(k_aH2P*y0/(1 + k_aH2P*y0)) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4c1(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn + k_txnBH*y2 + k_txnBH*y4*(k_aH2P*y0/(1 + k_aH2P*y0))
            - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*y2 + k_txnBH*y4*(k_aH2P*y0/(1 + k_aH2P*y0)) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

# =============================================================================
# MODEL 1A
# =============================================================================

def HBS_1a1A(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*y2 + k_txnBH*y4*k_aH2P*y0 - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4b1A(y, t, v):
    
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn + k_txnBH*y2 + k_txnBH*y4*k_aH2P*y0 - k_dR*y1
            - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*y2 + k_txnBH*y4*k_aH2P*y0 - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4c1A(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn + k_txnBH*y2 + k_txnBH*y4*k_aH2P*y0 - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*y2 + k_txnBH*y4*k_aH2P*y0 - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

# =============================================================================
# MODEL 1B
# =============================================================================

def HBS_1a1B(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*(y2 + y4) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4b1B(y, t, v):
    
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn + k_txnBH*(y2 + y4) - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*(y2 + y4) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4c1B(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y2 - k_dH1P*y0*y2,
            k_txn + k_txnBH*(y2 + y4) - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4,
            k_txnBH*(y2 + y4) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

# =============================================================================
# MODEL 1C
# =============================================================================

def HBS_1a1C(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y0*y4,
            k_txnBH*(y2 + y4) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4b1C(y, t, v):
    
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn + k_txnBH*(y2 + y4) - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y0*y2,
            k_txn - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y0*y4,
            k_txnBH*(y2 + y4) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

def HBS_4c1C(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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

    dydt = [k_basal + k_pM0*(y2 + y4) - k_dM0*y0,
            k_txn - k_dR*y1 - k_dH1R*y0*y1,
            k_tln*y1 - k_dP*y2 - k_dHP*O2_rate*y0*y2,
            k_txn + k_txnBH*(y2 + y4) - k_dR*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y0*y4,
            k_txnBH*(y2 + y4) - k_dR*y5,
            k_tln*y5 - k_dRep*y6]
        
    return dydt

# =============================================================================
# MODEL 2
# =============================================================================

def HBS_1a2(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn + k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) -
            k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn  + k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) -
            k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 2A
# =============================================================================

def HBS_1a2A(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*y3 + k_txnBH*y5*k_aH2P*y1 - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2A(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn + k_txnBH*y3 + k_txnBH*y5*k_aH2P*y1 - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*y5*k_aH2P*y1 - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2A(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn  + k_txnBH*y3 + k_txnBH*y5*k_aH2P*y1 -
            k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*y5*k_aH2P*y1 - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 2B
# =============================================================================

def HBS_1a2B(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2B(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn + k_txnBH*(y3 + y5) - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2B(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn  + k_txnBH*(y3 + y5) - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 2C
# =============================================================================

def HBS_1a2C(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y1*y5,
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2C(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn + k_txnBH*(y3 + y5) - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y1*y5,
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2C(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn  + k_txnBH*(y3 + y5) - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y1*y5,
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 2D
# =============================================================================

def HBS_1a2D(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y1*y5,
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2D(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn + k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1))
            - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y1*y5,
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2D(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn  + k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y1*y5,
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 2E
# =============================================================================

def HBS_1a2E(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5)) - 
            k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2E(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn + k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5))
            - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5)) - 
            k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2E(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txn  + k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5))
            - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5)) -
            k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 2F
# =============================================================================

def HBS_1a2F(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_pM1*y0 - k_dM1*y1 - k_aH2P*y1*y5,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5)) - 
            k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2F(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_pM1*y0 - k_dM1*y1 - k_aH2P*y1*y5,
            k_txn + k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5))
            - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5)) - 
            k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2F(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_pM1*y0 - k_dM1*y1 - k_aH2P*y1*y5,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn  + k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5))
            - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + y5)) -
            k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 2G
# =============================================================================

def HBS_1a2G(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_pM1*y0 - k_dM1*y1 - k_aH2P*y1*y5,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + k_pH2R*y5)) 
            - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b2G(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_pM1*y0 - k_dM1*y1 - k_aH2P*y1*y5,
            k_txn + k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + k_pH2R*y5))
            - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + k_pH2R*y5))
            - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4c2G(y, t, v):
 
    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_pM1*y0 - k_dM1*y1 - k_aH2P*y1*y5,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y1*y3,
            k_txn  + k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + k_pH2R*y5))
            - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y1*y5/(1 + k_aH2P*y1*y5 + k_pH2R*y5))
            - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

# =============================================================================
# MODEL 3
# =============================================================================

def HBS_1a3(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*y3 + k_txnBH*y5*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b3(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*y4 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dHP*O2_rate*y4 - k_dH1P*y1*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2_rate*y6,
            k_txnBH*y4 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
        
    return dydt

def HBS_4c3(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*y3 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dHP*O2_rate*y6,
            k_txnBH*y3 + k_txnBH*y6*(k_aH2P*y1/(1 + k_aH2P*y1)) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
    
    return dydt

# =============================================================================
# MODEL 3B
# =============================================================================

def HBS_1a3B(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b3B(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*(y4 + y6) - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dHP*O2_rate*y4 - k_dH1P*y1*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2_rate*y6,
            k_txnBH*(y4 + y6) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
        
    return dydt

def HBS_4c3B(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*(y3 + y6) - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dHP*O2_rate*y6,
            k_txnBH*(y3 + y6) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
    
    return dydt

# =============================================================================
# MODEL 3C
# =============================================================================

def HBS_1a3C(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y1*y5,
            k_txnBH*(y3 + y5) - k_dR*y6,
            k_tln*y6 - k_dRep*y7]
        
    return dydt

def HBS_4b3C(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_txnBH*(y4 + y6) - k_dR*y3,
            k_tln*y2 + k_tln2*y3 - k_dP*y4 - k_dHP*O2_rate*y1*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2_rate*y1*y6,
            k_txnBH*(y4 + y6) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
        
    return dydt

def HBS_4c3C(y, t, v):

    [[k_pM0, k_dM0, k_pM1, k_dM1, k_dH1R, k_dH1P, k_pH2R, k_dHP, k_aH2P, k_tln2, deg_ratio], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

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
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y1*y3,
            k_txn - k_dR*y4,
            k_txnBH*(y3 + y6) - k_dR*y5,
            k_tln*y4 + k_tln2*y5 - k_dP*y6 - k_dHP*O2_rate*y1*y6,
            k_txnBH*(y3 + y6) - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
    
    return dydt


# =============================================================================
# MODEL 4
# =============================================================================

def HBS_1a4(y, t, v):

    [[k_pH0, k_dH0, k_pH1, k_dH1, k_txnaH1, k_dH1R, k_dH1P, k_dHP, k_b, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    k_txnH = k_txnBH

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF, y1 = antisense HIF1a RNA,
    # y2 = HIF1a mRNA, y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = HIF2a protein,
    # y6 = active HIF2a protein, y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_basal + k_pH0*(y3 + y5) - k_dH0*y0 - k_b*y0*y5,
            k_txnH*y3 + k_txnH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5))
            - k_dR*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y0*y3,  
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5 - k_b*y0*y5,
            k_b*y0*y5 - k_dP*y6,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5)) 
            - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
        
    return dydt

def HBS_4b4(y, t, v):

    [[k_pH0, k_dH0, k_pH1, k_dH1, k_txnaH1, k_dH1R, k_dH1P, k_dHP, k_b, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    k_txnH = k_txnBH

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF, y1 = antisense HIF1a RNA,
    # y2 = HIF1a mRNA, y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = HIF2a protein,
    # y6 = active HIF2a protein, y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_basal + k_pH0*(y3 + y5) - k_dH0*y0 - k_b*y0*y5,
            k_txnH*y3 + k_txnH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5))
            - k_dR*y1,
            k_txn + k_txnBH*y3 + k_txnBH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5))
            - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y0*y3,  
            k_txn - k_dR*y4,
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5 - k_b*y0*y5,
            k_b*y0*y5 - k_dP*y6,
            
            k_txnBH*y3 + k_txnBH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5)) 
            - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
        
    return dydt

def HBS_4c4(y, t, v):

    [[k_pH0, k_dH0, k_pH1, k_dH1, k_txnaH1, k_dH1R, k_dH1P, k_dHP, k_b, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    k_txnH = k_txnBH

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF, y1 = antisense HIF1a RNA,
    # y2 = HIF1a mRNA, y3 = HIF1a protein, y4 = HIF2a mRNA, y5 = HIF2a protein,
    # y6 = active HIF2a protein, y7 = DsRED2 mRNA, y8 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8 = y

    dydt = [k_basal + k_pH0*(y3 + y5) - k_dH0*y0 - k_b*y0*y5,
            k_txnH*y3 + k_txnH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5))
            - k_dR*y1,
            k_txn - k_dR*y2 - k_dH1R*y1*y2,
            k_tln*y2 - k_dP*y3 - k_dHP*O2_rate*y3 - k_dH1P*y0*y3,
            k_txn + k_txnBH*y3 + k_txnBH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5))
            - k_dR*y4,
            
            k_tln*y4 - k_dP*y5 - k_dHP*O2_rate*y5 - k_b*y0*y5,
            k_b*y0*y5 - k_dP*y6,
            k_txnBH*y3 + k_txnBH*(k_aH2P*y6/(1 + k_aH2P*y6 + k_iH2P*y5)) 
            - k_dR*y7,
            k_tln*y7 - k_dRep*y8]
        
    return dydt

# =============================================================================
# MODEL 4A
# =============================================================================

def HBS_1a4A(y, t, v):

    [[k_pH0, k_dH0, k_pH1, k_dH1, k_txnaH1, k_dH1R, k_dH1P, k_dHP, k_b, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    k_txnH = k_txnBH

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF, y1 = active HAF,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = active HIF2a protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_basal + k_pH0*(y4 + y6) - k_dH0*y0 - k_pH1*y0,
            k_pH1*y0 - k_dH1*y1 - k_b*y1*y6,
            k_txnH*y4 + k_txnH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6))
            - k_dR*y2,
            k_txn - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4 - k_dH1P*y1*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2_rate*y6 - k_b*y1*y6,
            k_b*y1*y6 - k_dP*y7,
            k_txnBH*y4 + k_txnBH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6)) 
            - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

def HBS_4b4A(y, t, v):

    [[k_pH0, k_dH0, k_pH1, k_dH1, k_txnaH1, k_dH1R, k_dH1P, k_dHP, k_b, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    k_txnH = k_txnBH

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF, y1 = active HAF,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = active HIF2a protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_basal + k_pH0*(y4 + y6) - k_dH0*y0 - k_pH1*y0,
            k_pH1*y0 - k_dH1*y1 - k_b*y1*y6,
            k_txnH*y4 + k_txnH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6))
            - k_dR*y2,
            k_txn + k_txnBH*y4 + k_txnBH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6))
            - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4 - k_dH1P*y1*y4,
            k_txn - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2_rate*y6 - k_b*y1*y6,
            k_b*y1*y6 - k_dP*y7,
            k_txnBH*y4 + k_txnBH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6)) 
            - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

def HBS_4c4A(y, t, v):

    [[k_pH0, k_dH0, k_pH1, k_dH1, k_txnaH1, k_dH1R, k_dH1P, k_dHP, k_b, k_aH2P, k_iH2P], O2] = v

    #parameters that will be held constant:
    k_txnBH = 1.0
    k_dR = 2.7 #1/h
    k_tln = 1 #1/h
    k_dP = 0.35 #1/h
    k_dRep = 0.029 #1/hr

    k_basal = 1.0 #au (MARS with NEW Mechanism)
    k_txn = 1 #au/h
    k_txnH = k_txnBH

    #when in normoxic incubator, use normoxic pO2
    if O2 == 138:
        O2_rate = O2

    #Exponential decrease in pO2 from fit to diffusion equation
    else:
        # O2_rate = max(144.41*np.exp(-0.011*(60*t)), O2)
        O2_rate = (135.81 - 7.6)*np.exp(-(3600*t-500)*4.25e-4) + 7.6
        if O2_rate > 138:
            O2_rate = 138
        
    # y holds these state variables: y0 = HAF, y1 = active HAF,
    # y2 = antisense HIF1a RNA, y3 = HIF1a mRNA, y4 = HIF1a protein,
    # y5 = HIF2a mRNA, y6 = HIF2a protein, y7 = active HIF2a protein,
    # y8 = DsRED2 mRNA, y9 = DsRED2 protein

    y0, y1, y2, y3, y4, y5, y6, y7, y8, y9 = y

    dydt = [k_basal + k_pH0*(y4 + y6) - k_dH0*y0 - k_pH1*y0,
            k_pH1*y0 - k_dH1*y1 - k_b*y1*y6,
            k_txnH*y4 + k_txnH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6))
            - k_dR*y2,
            k_txn - k_dR*y3 - k_dH1R*y2*y3,
            k_tln*y3 - k_dP*y4 - k_dHP*O2_rate*y4 - k_dH1P*y1*y4,
            k_txn  + k_txnBH*y4 + k_txnBH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6))
            - k_dR*y5,
            k_tln*y5 - k_dP*y6 - k_dHP*O2_rate*y6 - k_b*y1*y6,
            k_b*y1*y6 - k_dP*y7,
            k_txnBH*y4 + k_txnBH*(k_aH2P*y7/(1 + k_aH2P*y7 + k_iH2P*y6)) 
            - k_dR*y8,
            k_tln*y8 - k_dRep*y9]
        
    return dydt

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
