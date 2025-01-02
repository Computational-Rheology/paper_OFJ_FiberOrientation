#%% Modules
import numpy as np
from scipy.integrate import solve_ivp

from orientation_models import create_model
from closure_models import create_closure

import pandas as pd
from inspect import signature
import os

#%% Velocity gradient func
def vel_grad(t):
        return np.array([[0, 1.0, 0],
                         [0, 0, 0],
                         [0, 0, 0]])
    
#%% Check function signature and keyword arguments
def parse_func_args(func, kwargs):
        func_params = str(signature(func)).replace("(t, A2,", "").replace(")","").replace(" ", "").split(",")
        params = []
        for p in func_params:
            if p in kwargs.keys():               
                params.append(kwargs[p])
                kwargs.pop(p, None)
        
        if len(params) != len(func_params):
            raise RuntimeError(f"problem with the parameters: {kwargs.keys()}")
            
        return params

#%% Create data folder
def create_folder(name):
    abs_path = os.path.abspath(name)
    if not os.path.exists(abs_path):
        os.mkdir(abs_path)
    return True

#%%  Main
def main(t, 
         A2_0, 
         models_and_kwargs, 
         closure_list, 
         save_data=False,
         save_folder=None):
       
    for name, kwargs in models_and_kwargs.items():
        model_fn = create_model(name)
        
        for closure_name in closure_list:
            kwargs["closure_fn"] = create_closure(closure_name)
            args_to_func = parse_func_args(model_fn, kwargs.copy())
            
            sol = solve_ivp(    model_fn,
                                (t[0], t[-1]), 
                                A2_0.ravel(),
                                method = 'RK45',
                                t_eval = t,
                                args=(*args_to_func,),
                                rtol = 1e-12,
                                atol = 1e-12
                            )
            
            df = pd.DataFrame({
                               'time': t,
                               'A2_0': sol.y[0, :], # xx
                               'A2_1': sol.y[1, :], # xy
                               'A2_2': sol.y[2, :], # xz
                               'A2_3': sol.y[4, :], # yy
                               'A2_4': sol.y[5, :], # yz
                               'A2_5': sol.y[8, :], # zz
                              })
        
            if save_data:
                df.to_csv(f'{save_folder}/{name}_{closure_name}.csv', index=False)
            
#%% Run
if __name__ == '__main__':
    
    save_folder = "./data_check"
    
    create_folder(save_folder)

    A2_0 = np.eye(3)*(1.0/3) 
    
    dt = 0.01
    
    t = np.arange(0, 65.0 + dt, dt)
    
    models_and_kwargs = {
                            "FT": {
                                    "vel_grad": vel_grad, 
                                    "CI": 0.0311, 
                                    "xi": 1.0
                                  },
                            
                            "iARD_RPR": {
                                            "vel_grad": vel_grad, 
                                            "CI": 0.0562, 
                                            "xi": 1.0, 
                                            "CM": 0.9977, 
                                            "alpha": 0.0, 
                                            "beta": 0.0
                                        },
                            
                            "MRD": {
                                        "vel_grad": vel_grad, 
                                        "CI": 0.0198, 
                                        "xi": 1.0, 
                                        "D1": 1.0, 
                                        "D2": 0.4796,
                                        "D3": 0.0120
                                    }
                         }
    
    closures=['hybrid', 'IBOF', 'ORE']
    
    main(t, A2_0, models_and_kwargs, closures, True, save_folder)    
    
    