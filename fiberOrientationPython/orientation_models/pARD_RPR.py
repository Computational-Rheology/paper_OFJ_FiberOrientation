#%% Mods
import numpy as np
from .flow_tensors import flow_tensors

#%% Funcs
def pARD_RPR(t, A2, vel_grad, closure_fn, xi, CI, Omega, alpha):
    """ 
        principal Anisotropic Rotary Diffusion model (pARD-RPR)
        
        Parameters
        ----------
            t : scalar
                time
                
            A2 : (9,) numpy array
                1D fiber orientation tensor.
                
            vel_grad : function 
                Velocity gradient
        
            closure_fn: function 
                Function to compute A4 given A2
            
            xi : scalar
                Shape factor
        
            CI : scalar
                fiber-fiber interaction coefficient
        
            Omega: scalar [should be between 0.5 and 1]
        
            alpha : scalar [between 0 and 1]
                strong factor of the transient orientation slow-down rate
    
        Literature:
            H.-C. Tseng, R.-Y. Chang, and C.-H. Hsu, 
            “The use of principal spatial tensor to predict anisotropic fiber orientation in concentrated fiber suspensions,”
            J. Rheol. (N. Y. N. Y)., vol. 62, no. 1, pp. 313–320, 2018, doi: 10.1122/1.4998520.
    """
    
    A2 = A2.reshape(3, 3)
    A4 = closure_fn(A2)    
    D, W, shrRate = flow_tensors(t, vel_grad)
        
    ## Calculates Adot_HD -> Jeffery hydrodynamics
    Adot_HR =(
                (np.tensordot(W, A2, axes=1) - np.tensordot(A2, W, axes=1))
               + xi*(
                        np.tensordot(D, A2, axes=1) + np.tensordot(A2, D, axes=1)
                        - 2.0*np.tensordot(A4, D, axes=2)
                    )
             )
    

    lambda_, R_ = np.linalg.eig(A2)
    idx = lambda_.argsort()[::-1]   
    lambda_ = lambda_[idx]
    R_ = R_[:,idx]
    
    # Defines C_hat
    C_hat = np.array([[1.0, 0 , 0],
                      [0, Omega, 0 ],
                      [0, 0, 1.0 - Omega]
                      ])
    # Calculate C
    C = CI*(np.linalg.multi_dot([R_, C_hat, R_.T]))
    
    # # or
    # C=CI*( np.tensordot(R_[:,0],R_[:,0],axes=0)
    #       +Omega*np.tensordot(R_[:,1],R_[:,1],axes=0)
    #       +(1-Omega)*np.tensordot(R_[:,2],R_[:,2],axes=0)
    #       )
    
    # Calculate Adot_pARD
    Adot_pARD = shrRate*(
                           2.0*C
                          -2.0*np.trace(C)*A2
                          -5.0*np.tensordot(C, A2, axes=1)
                          -5.0*np.tensordot(A2, C, axes=1)
                          +10.0*np.tensordot(A4, C,axes=2)
                        )
    
    # Calculate Adot_RPR
    Adot_HR_pARD= Adot_HR + Adot_pARD
    
    # Calculate lambda_dot through IOK (Intrinsic Orientation Kinetic mechanism)
    lambda_Adot_HR_pARD = np.diag(np.linalg.multi_dot([R_.T, Adot_HR_pARD, R_]))
    
    LambdaDot_IOK = np.zeros([3,3])
    for i in range(3):
        LambdaDot_IOK[i, i] = alpha*(lambda_Adot_HR_pARD[i%3])
                            
    
    # Calculate Adot_RPR
    Adot_RPR = np.linalg.multi_dot([-R_, LambdaDot_IOK, R_.T])
     
    result = (Adot_HR_pARD + Adot_RPR)
    
    return result.ravel()
 