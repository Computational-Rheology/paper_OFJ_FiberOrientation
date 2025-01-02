#%% Mods
import numpy as np
from .flow_tensors import flow_tensors
        
#%% Func  
def MRD(t, A2, vel_grad, closure_fn, xi, CI, D1, D2, D3):
    """ 
        Moldflow rotary diffusion model
        
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
        
            D1 : scalar
                coefficient. Show how the rotational diffusion is biased 
                towards the direction of the principal vectors of fiber orientation tensor  
                    
            D2 : scalar
                coefficient. Show how the rotational diffusion is biased 
                towards the direction of the principal vectors of fiber orientation tensor    
                            
            D3 : scalar
                coefficient. Show how the rotational diffusion is biased 
                towards the direction of the principal vectors of fiber orientation tensor   
            
        
        
        Literature:
            A. Bakharev, H. Yu, S. Ray, R. Speight, and J. Wang, 
            “Using new anisotropic rotational diffusion model to improve prediction of short fibers in thermoplastic injection molding,”
            Annu. Tech. Conf. - ANTEC, Conf. Proc., vol. 2018-May, 2018.
            
            A. J. Favaloro and C. L. Tucker, 
            “Analysis of anisotropic rotary diffusion models for fiber orientation,”
            Compos. Part A Appl. Sci. Manuf., vol. 126, no. August, 2019, doi: 10.1016/j.compositesa.2019.105605. 
            *Note: Uses a different formulation for this model (Phelps-Tucker formulation)
            
           S. K. Kugler, G. M. Lambert, C. Cruz, A. Kech, T. A. Osswald, and D. G. Baird,
           “Macroscopic fiber orientation model evaluation for concentrated short fiber reinforced polymers in comparison to experimental data,” 
           Polym. Compos., vol. 41, no. 7, pp. 2542–2556, 2020, doi: 10.1002/pc.25553.
    """

    # The fiber orientation model is defined as:
    # dA/dt =  Adot_HD + Adot_MRD;

    A2 = A2.reshape(3, 3)
    A4 = closure_fn(A2)
    
    D, W, shrRate = flow_tensors(t, vel_grad)

    # Calculate Adot_HD (Jeffery Hydrodynamics)
    Adot_HD = (
                 (np.tensordot(W, A2, axes=1) - np.tensordot(A2, W, axes=1))
                 + xi*(
                         np.tensordot(D, A2, axes=1) + np.tensordot(A2, D, axes=1)
                         - 2.0*np.tensordot(A4, D, axes=2)
                      )
              )
    
    # eigendecomposition of  A2 
    lambda_, e_ = np.linalg.eig(A2)
    
    # Sort eigenvalues in descending order. Sort eigenVectors in order of the eigenValues
    idx = lambda_.argsort()[::-1]   
    lambda_ = lambda_[idx]
    e_ = e_[:, idx]

    # C_MRD = CI*( D1*np.tensordot(e_[:,0],e_[:,0], axes=0)
    #             +D2*np.tensordot(e_[:,1],e_[:,1], axes=0)
    #             +D3*np.tensordot(e_[:,2],e_[:,2], axes=0)
    #             )
    
    D = np.array([
                    [D1, 0, 0],
                    [0, D2, 0],
                    [0, 0, D3]
                  ])
    
    C_MRD = CI*(np.linalg.multi_dot([e_, D, e_.T]))   
     
    ## Original Formulation
    # Adot_MRD = 2.0*shrRate*(C_MRD - np.trace(C_MRD)*A2)
    
    # From Favaloro and Tucker
    Adot_MRD = shrRate*( 2.0*C_MRD
                        -2.0*np.trace(C_MRD)*A2
                        -5.0*(np.tensordot(C_MRD,A2, axes=1) + np.tensordot(A2, C_MRD, axes=1))
                        +10.0*np.tensordot(A4, C_MRD, axes=2))
        
    result = Adot_HD + Adot_MRD
    

    return result.ravel()
