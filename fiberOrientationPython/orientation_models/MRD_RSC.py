#%% Modules
import numpy as np
from .flow_tensors import flow_tensors

#%% Funcs          
def MRD_RSC(t, A2, vel_grad, closure_fn, xi, CI, k, D1, D2, D3):
    """ 
         Moldflow rotary diffusion model with Reduced Strain Closure (MRD-RSC)
        
        Parameters
        ----------
            t : scalar
                time
                
            A2 : (9,) numpy array
                1D fiber orientation tensor.
                
            vel_grad : function 
                Function to compute velocity gradient at time t
        
            closure_fn: function 
                Function to compute A4 given A2
            
            xi : scalar
                Shape factor
        
            CI : scalar
                fiber-fiber interaction coefficient    
                
            k : scalar
                 orientation kinetic slow-down parameter   
        
            D1 : scalar
                coefficient of asymmetry. Show how the rotational diffusion is biased 
                towards the direction of the principal vectors of fiber orientation tensor  
                    
            D2 : scalar
                coefficient of asymmetry. Show how the rotational diffusion is biased 
                towards the direction of the principal vectors of fiber orientation tensor    
                            
            D3 : scalar
                coefficient of asymmetry. Show how the rotational diffusion is biased 
                towards the direction of the principal vectors of fiber orientation tensor   
            
        Literature:
            S. K. Kugler, G. M. Lambert, C. Cruz, A. Kech, T. A. Osswald, and D. G. Baird, 
            “Macroscopic fiber orientation model evaluation for concentrated short fiber reinforced 
            polymers in comparison to experimental data,” 
            Polym. Compos., vol. 41, no. 7, pp. 2542–2556, 2020, doi: 10.1002/pc.25553.
    """

    D, W, shrRate = flow_tensors(t, vel_grad)
    A2 = A2.reshape(3, 3)
    A4 = closure_fn(A2)
    
    # spectral decomposition of A2
    lambda_, e_ = np.linalg.eig(A2)
    
    # Sort eigenvalues in descending order.
    idx = lambda_.argsort()[::-1]   
    lambda_ = lambda_[idx]
    e_ = e_[:,idx]

    C_hat = np.array([
                        [D1, 0, 0],
                        [0, D2, 0],
                        [0, 0, D3]
                     ])

    C_MRD = CI*np.linalg.multi_dot([e_, C_hat, e_.T])

    # Calculate M and L
    M= np.zeros([3,3,3,3])
    L= np.zeros([3,3,3,3])
    
    for i in np.arange(3):
        M += np.tensordot(np.tensordot(e_[:, i], e_[:, i], axes=0),
                          np.tensordot(e_[:, i], e_[:, i], axes=0), axes=0)
        
        L += lambda_[i]*np.tensordot(np.tensordot(e_[:, i], e_[:, i], axes=0),
                                     np.tensordot(e_[:, i], e_[:, i], axes=0), axes=0)
   
    result = (
                 (np.tensordot(W, A2, axes=1) - np.tensordot(A2, W, axes=1))
                +xi*(
                       np.tensordot(D, A2, axes=1) + np.tensordot(A2, D, axes=1)
                      -2.0*np.tensordot((A4 + (1 - k)*(L - np.tensordot(M, A4, axes=2))), D, axes=2)
                    
                    )
                + 2.0*shrRate*k*(C_MRD - np.trace(C_MRD)*A2)
             )
    
    
    return result.ravel()
