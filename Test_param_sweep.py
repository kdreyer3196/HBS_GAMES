# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 18:41:49 2022

@author: Katie_Dreyer
"""

p_ref = [3.5, 0.24]

real_param_labels_all = ['p1', 'p2']

orders_of_mag = [-3, -2, -1, 0, 1, 2, 3]

param_sweeps = {}

for index, param in enumerate(p_ref):
    param_label = real_param_labels_all[index]
    
    param_sweeps[param_label] = []
    
    
    for exp in orders_of_mag:
        param_new = param*10**exp
        
        param_sweeps[param_label].append(param_new)  
    
# print(param_sweeps)


param_lists = {}

for index, param in enumerate(p_ref):
    
    param_label = real_param_labels_all[index]
    
    param_lists[param_label + ' sweep'] = []
        
    for p_new in param_sweeps[param_label]:
        
        new_params = p_ref.copy()

        
        new_params[index] = p_new
        
        # print(new_params)
        
        param_lists[param_label + ' sweep'].append(new_params)
                
        
print(param_lists)
    

