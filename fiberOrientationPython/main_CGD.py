#%% Modules
import numpy as np
from scipy.integrate import solve_ivp
from closure_models import create_closure

import pandas as pd
from tqdm import tqdm

#%% funcs
def FT_r(r, A2, z, closure_fn, xi, CI, b, Q):
    A2 = A2.reshape(3, 3)

    A4 = closure_fn(A2)
    
    U = vel(r, z, b, Q)
    
    L = vel_grad(r, z, b, Q)
    
    D = 0.5*(L + L.T)
    
    W = 0.5*(L - L.T)
        
    shrRate = np.sqrt(2.0*np.tensordot(D, D, axes=2))
       
    I = np.eye(3)
    
    result =   ( 
                  ((W @ A2) - (A2 @ W))
                  +xi*(
                          (D @ A2) + (A2 @ D)
                          - 2.0*np.tensordot(A4, D, axes=2)
                      )
                  + 2.0*(CI*shrRate)*(I - 3.0*A2) 
               )/U
    
    return result.ravel()


def vel(r, z, b, Q): 
    aux1 = (3.0*Q)/(8.0*np.pi*r*b)
    aux2 = (1.0-(z**2/b**2))
    
    return aux1*aux2


def vel_grad(r, z, b, Q):
    aux = (3.0*Q)/(8.0*np.pi*r*b)
    xx = (1.0/r)*(1.0-(z**2/b**2))
    yy = (1.0/r)*(1.0-(z**2/b**2))
    xz = (2.0/b)*(z/b)
    
    return aux*np.array([[-xx,  0.0,   -xz],
                         [0.0,  yy,  0.0],
                         [0.0, 0.0, 0.0]])

#%%  Main
def main():

    closure_name = "IBOF"
    
    closure_fn = create_closure(closure_name)
    
    xi = 1.0
    
    CI = 0.001
    
    b = 1.5/1000
    
    Q = 0.000134774
    
    dz = 2.34375e-05
    
    z = np.arange(0.00001171875, 0.00298828125 + dz, dz) 
    
    z = z - b

    r = np.linspace(0.01, 0.12, 3000)
        
    A2_0 = np.eye(3)*(1.0/3)  # Initial A2 in isotropic state 
    
    df = pd.DataFrame(np.zeros((len(r)*len(z), 8)))
    
    df.columns=['r','z',
                'A_{11}', 'A_{12}', 'A_{13}', 'A_{22}', 'A_{23}', 'A_{33}']
    
    # Populate r and z
    for i in range(len(z)):
        df.iloc[i*len(r):(i+1)*len(r), 0] = r
        df.iloc[i*len(r):(i+1)*len(r), 1] = z[i]
    
    # Solve ODE
    # Loops through every thickness and integrates over the radius
    
    for i in tqdm(range(len(z))): 
        sol = solve_ivp(
                            FT_r, 
                            (r[0], r[-1]), 
                            A2_0.ravel(),
                            method='RK45', 
                            t_eval = r, 
                            args=(z[i], closure_fn, xi, CI, b, Q),
                            atol=1e-12,
                            rtol=1e-12
                        )
        
        df.iloc[i*len(r):(i + 1)*len(r), 2:] = sol.y.T[:, [0, 1, 2, 4, 5, 8]]
    
    df.to_csv(closure_name + '_' + 'disk' + '.csv', index=False)     # Save data to csv
     
    return df
    
#%% Run
if __name__ == '__main__':
    df = main()
    