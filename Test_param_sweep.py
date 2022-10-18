# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 18:41:49 2022

@author: Katie_Dreyer
"""
def get_param_sweeps(p_ref, orders_of_mag):
    
    param_sweeps = {}
    
    for index, param in enumerate(p_ref):
        param_label = real_param_labels_all[index]
        
        param_sweeps[param_label] = []
        
        
        for exp in orders_of_mag:
            param_new = param*10**exp
            
            param_sweeps[param_label].append(param_new)  
        
    return param_sweeps

def get_param_lists(p_ref, param_sweeps):
    
    param_lists = {}
    
    for index, param in enumerate(p_ref):
        
        param_label = real_param_labels_all[index]
        
        param_lists[param_label + ' sweep'] = []
            
        for p_new in param_sweeps[param_label]:
            
            new_params = p_ref.copy()
    
            new_params[index] = p_new
                        
            param_lists[param_label + ' sweep'].append(new_params)
            
    return param_lists





orders_of_mag = [-3, -2, -1, 0, 1, 2, 3]

p_ref = [3.5, 0.24]

real_param_labels_all = ['p1', 'p2']

param_sweeps = get_param_sweeps(p_ref, orders_of_mag)

param_lists = get_param_lists(p_ref, param_sweeps)

# print(param_sweeps, param_lists)
    
for param_label, param_sets in param_lists.items():
    print(param_label, ': ', param_sets)
    for params in param_sets:
        print(params)



#k_p = 0.02

# t = np.linspace(0,120,31)

# k_p = np.piecewise(t, [t < 44, t >= 44], [0, k_p])

# print(k_p)