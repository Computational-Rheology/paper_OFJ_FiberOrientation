#%% Mods
import numpy as np

#%% Funcs
def flow_tensors(t, vel_grad):
    """    
        Parameters
        ----------
        t : scalar
            time
        vel_grad : function
            Function to compute the velocity gradient.
            Must return a 3x3 Numpy array. The gradient should be defined as 
            \nabla \mathbf{v}_{ij} = \frac{\partial v_{i}}{\partial x_{j}}
    
        Returns
        -------
        D : 3x3 Numpy array
            Rate-of-deformation tensor [symmetric part of velocity gradient]
        W : 3x3 Numpy array
            Vorticity tensor [anti-symmetric part of velocity gradient]
        shrRate : scalar
            scalar magnitude of the velocity gradient (shear-rate)

    """
    D = 0.5*(vel_grad(t) + vel_grad(t).T)
    
    W = 0.5*(vel_grad(t) - vel_grad(t).T)
        
    shrRate = np.sqrt(2*np.tensordot(D, D, axes=2))
 
    return (D, W, shrRate)