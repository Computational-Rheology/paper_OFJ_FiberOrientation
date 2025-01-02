#%% Mods
import numpy as np
from .flow_tensors import flow_tensors

#%% Func
def iARD_RPR(t, A2, vel_grad, closure_fn, xi, CI, CM, alpha, beta):
    
    """ 
        Objective iARD-RPR Model
        
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
                form factor
        
            CI : scalar
                fiber-fiber interaction coefficient
        
            CM : scalar
                fiber-matrix interaction coefficient
                
            alpha : scalar
                strong factor of slow-down rate
                
            beta : scalar 
                weak factor of slow-down rate
            
        Literature:
            H.-C. Tseng, R.-Y. Chang, and C.-H. Hsu, 
            “An objective tensor to predict anisotropic fiber orientation in concentrated suspensions,”
            J. Rheol. (N. Y. N. Y)., vol. 60, no. 2, pp. 215–224, 2016, 
            doi: https://doi.org/10.1122/1.4939098.
        
    """

    # In  this model the fiber orientation model is defined as:
    # dA/dt =  Adot_HD + Adot_iARD(CI, CM) + Adot_RPR(alpha)

    A2 = A2.reshape(3, 3)
    A4 = closure_fn(A2)
    
    D, W, shrRate = flow_tensors(t, vel_grad)
     
    ## Calculates Adot_HD [Jeffery Hydrodynamics]
    Adot_HR = (
                    (np.tensordot(W, A2, axes=1) - np.tensordot(A2, W, axes=1))
                    + xi*(
                           np.tensordot(D, A2, axes=1) + np.tensordot(A2, D, axes=1)
                           - 2.0*np.tensordot(A4, D, axes=2)
                         )
               )
    
    ## Calculate Adot_iARD
    D_square = np.tensordot(D, D, axes=1)
    
    # Calculate the norm of D²
    D_square_norm = np.sqrt(0.5*np.tensordot(D_square, D_square, axes=2))
    
    # Calculate the fiber-rotary-resistance
    L_tilda = D_square/D_square_norm
    
    # Create the 2nd order identity tensor   
    I = np.eye(3)
    
    # Calculate the ARD tensor
    D_r= CI*(I - CM*L_tilda)
    
    # Calculate Adot_iARD
    Adot_iARD = shrRate*(
                           2.0*D_r
                          -2.0*np.trace(D_r)*A2
                          -5.0*np.tensordot(D_r, A2, axes=1)
                          -5.0*np.tensordot(A2, D_r, axes=1)
                          +10.0*np.tensordot(A4, D_r, axes=2)
                        )
       
    ## Calculate Adot_RPR
    
    # Groups the HD and iARD components 
    Adot_HR_iARD= Adot_HR + Adot_iARD
    
    # Calculates eigenvalues and eigenvectors of A2
    # Note: possibly using "eigh" would be a better approach
    lambda_, e_ = np.linalg.eig(A2)
    
    # Sort eigenvalues in descending order. 
    idx = lambda_.argsort()[::-1]   
    lambda_ = lambda_[idx]
    e_ = e_[:, idx]

    # Calculates eigenvalues of the combined Adot_HR and Adot_iARD with the eigenvectors of A2
    lambda_Adot_HR_iARD = np.diag(np.linalg.multi_dot([e_.T, Adot_HR_iARD, e_]))
    
    # Calculates the eigenvalues according to the IOK (Intrinsic Orientation Kinetic mechanism)
    LambdaDot_IOK = np.zeros([3, 3])
    
    for i in range(3):
        LambdaDot_IOK[i,i] = alpha*(
                                     lambda_Adot_HR_iARD[i%3]
                                     - beta*(
                                              lambda_Adot_HR_iARD[i%3]**2
                                              + 2.0*lambda_Adot_HR_iARD[(i + 1)%3]*lambda_Adot_HR_iARD[(i + 2)%3]
                                            )
                            
                                    )
    
    # Calculates Adot_RPR
    Adot_RPR = np.linalg.multi_dot([-e_, LambdaDot_IOK, e_.T])
     
    # iARD-RPR objective fiber orientation equation
    result = (Adot_HR_iARD + Adot_RPR)
    
    return result.ravel()
