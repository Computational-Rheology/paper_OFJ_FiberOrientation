#%% Modules
import numpy as np
from .flow_tensors import flow_tensors

#%% Fun
def FT_RSC(t, A2, vel_grad, closure_fn, xi, CI, k):
    """
        Folgar-Tucker model with Reduced-Strain Closure Model
        
       Parameters
        ----------
            t : scalar
                time
                
            A2 : (9,) numpy array
                1D fiber orientation tensor
                
            vel_grad : function 
                Velocity gradient
        
            closure_fn: function 
                Function to compute A4 given A2
            
            xi : scalar
                Shape factor
        
            CI : scalar
                fiber-fiber interaction coefficient
        
            k : scalar
                Slow-down parameter
        
        Literature:
                J. Wang, J. F. O’Gara, and C. L. Tucker, 
                “An objective model for slow orientation kinetics in concentrated fiber suspensions: Theory and rheological evidence,” 
                J. Rheol. (N. Y. N. Y)., vol. 52, no. 5, pp. 1179–1200, 2008, doi: 10.1122/1.2946437
                 
                J. H. Phelps and C. L. Tucker, 
                “An anisotropic rotary diffusion model for fiber orientation in short- and long-fiber thermoplastics,”
                J. Nonnewton. Fluid Mech., vol. 156, no. 3, pp. 165–176, 2009, doi: 10.1016/j.jnnfm.2008.08.002.
    """

    A2 = A2.reshape(3, 3)
    A4 = closure_fn(A2)
    
    # Compute eigenvalues and eigenvectors of A
    lambda_, e_ = np.linalg.eig(A2)
    
    # Sorts eigenvalues in descending order
    idx = lambda_.argsort()[::-1]   
    lambda_ = lambda_[idx]
    e_ = e_[:,idx]
    
    # Calculate M and L
    M= np.zeros([3,3,3,3])
    L= np.zeros([3,3,3,3])
    
    for i in np.arange(3):
        M += np.tensordot(np.tensordot(e_[:,i],e_[:,i],axes=0), 
                          np.tensordot(e_[:,i],e_[:,i],axes=0), axes=0)
        
        L += lambda_[i]*np.tensordot(np.tensordot(e_[:,i],e_[:,i],axes=0), 
                                     np.tensordot(e_[:,i],e_[:,i],axes=0), axes=0)


    #     w, v = np.linalg.eig(a)
    # L = (
    #     w[0] * np.einsum("i,j,k,l->ijkl", v[:, 0], v[:, 0], v[:, 0], v[:, 0])
    #     + w[1] * np.einsum("i,j,k,l->ijkl", v[:, 1], v[:, 1], v[:, 1], v[:, 1])
    #     + w[2] * np.einsum("i,j,k,l->ijkl", v[:, 2], v[:, 2], v[:, 2], v[:, 2])
    # )
    # M = (
    #     np.einsum("i,j,k,l->ijkl", v[:, 0], v[:, 0], v[:, 0], v[:, 0])
    #     + np.einsum("i,j,k,l->ijkl", v[:, 1], v[:, 1], v[:, 1], v[:, 1])
    #     + np.einsum("i,j,k,l->ijkl", v[:, 2], v[:, 2], v[:, 2], v[:, 2])
    # )

    # tensor4 = A + (1.0 - kappa) * (L - np.einsum("ijmn,mnkl->ijkl", M, A))    
   
    # Calculates the rate-of-deformation, vorticity and shear-rate
    D, W, shrRate = flow_tensors(t, vel_grad)
    
    # Creates the 2nd order identity tensor   
    I = np.eye(3)
    
    # Combines 4th order tensor operations
    tmp1 = A4 + (1.0 - k)*(L - np.tensordot(M, A4, axes=2))
    
    # RSC fiber orientation equation
    result =(
              (np.tensordot(W, A2, axes=1) - np.tensordot(A2, W, axes=1))
              +xi*(
                     np.tensordot(D, A2, axes=1) + np.tensordot(A2, D, axes=1)
                     - 2.0*np.tensordot(tmp1, D, axes=2)
                  )
              + 2.0*k*CI*shrRate*(I - 3.0*A2) 
            )
    
    return result.ravel()
