#%% Modules
import numpy as np
from .flow_tensors import flow_tensors

#%% Funcs
def FT(t, A2, vel_grad, closure_fn, xi, CI):
    
    """ 
        Folgar-Tucker Model
        
        Parameters:
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
     
        Literature:
            Advani, S. G., & Tucker, C. L. (1987). 
            "The Use of Tensors to Describe and Predict Fiber Orientation in Short Fiber Composites." 
            Journal of Rheology, 31(8), 751–784. https://doi.org/10.1122/1.549945
            
            S. G. Advani and C. L. Tucker, 
            “Closure approximations for three‐dimensional structure tensors,”
            J. Rheol. (N. Y. N. Y)., vol. 34, no. 3, pp. 367–386, 1990, doi: 10.1122/1.550133.
            
            Favaloro, A. J., & Tucker, C. L. (2019). 
            "Analysis of anisotropic rotary diffusion models for fiber orientation."
            Composites Part A: Applied Science and Manufacturing, 126(August). https://doi.org/10.1016/j.compositesa.2019.105605
            
            Kugler, S. K., Kech, A., Cruz, C., & Osswald, T. (2020). 
            "Fiber Orientation Predictions—A Review of Existing Models."
            Journal of Composites Science, 4(2), 69. https://doi.org/10.3390/jcs4020069
    """
       
    # A2 enters the function in ravel format (flattened orientation tensor)
    # Reshape to 3x3
    A2 = A2.reshape(3, 3)
          
    # Computes A4 given A2 through a closure relationship
    A4 = closure_fn(A2)
    
    #%% Method 1: As proposed in Advani and Tucker (1987)

    # # Computes rate-of-deformation tensor from the velocity gradient
    # D = (np.transpose(vel_grad(t)) + vel_grad(t))
    
    # # Computes the vorticity tensor from the velocity gradient
    # W = (np.transpose(vel_grad(t)) - vel_grad(t))

    # # calculates the shear rate
    # shrRate=np.sqrt(0.5*np.tensordot(D, D, axes=2))
       
    # I = np.eye(3)

    # result =   (        
    # -(1.0/2)*(np.einsum("ik,kj->ij", W, A2)
    # - np.einsum("ik,kj->ij", A2, W))
    # + (1.0/2)*xi*(
    #     np.einsum("ik,kj->ij", D, A2)
    #     + np.einsum("ik,kj->ij", A2, D)
    #     - 2.0*np.einsum("kl,ijkl->ij", D, A4) )
    # + 2.0*(CI*shrRate*(I - 3.0*A2))
    # )
    # return result.ravel()

    #%% Method 2: As reported in Favaloro and Tucker (2019)
    
    # D, W, shrRate = flowFieldTensors(t, vel_grad)
       
    # I = np.eye(3)
    
    # result =   ( 
    #               (np.tensordot(W, A2, axes=1) - np.tensordot(A2, W, axes=1))
    #               +xi*(np.tensordot(D, A2, axes=1) + np.tensordot(A2, D, axes=1) 
    #                    - 2.0*np.tensordot(A4, D, axes=2))
    #               +6.0*(CI*shrRate)*(I/3 - A2) 
    #            )
    
    # return result.ravel()
    
    #%%    Method 3: As reported in Kugler et al (2020)     
    D, W, shrRate = flow_tensors(t, vel_grad)
       
    I = np.eye(3)
    
    result =   ( 
        np.einsum("ik,kj->ij", W, A2)
        - np.einsum("ik,kj->ij",A2, W)
        + xi
        * (
              np.einsum("ik,kj->ij", D, A2)
            + np.einsum("ik,kj->ij", A2, D)
            - 2.0*np.einsum("ijkl,kl->ij", A4, D)
        )
        + 2.0*CI*shrRate*(I - 3.0*A2)
    )

    return result.ravel()

