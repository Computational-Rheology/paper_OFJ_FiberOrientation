#%% Modules
import numpy as np

#%% Funcs
def linear(A2):
    I = np.eye(3)
    
    A4_linear =((1.0/7.0)*(  np.einsum("ij,kl", A2, I)
                           + np.einsum("ik,jl", A2, I)
                           + np.einsum("il,jk", A2, I)
                           + np.einsum("kl,ij", A2, I)
                           + np.einsum("jl,ik", A2, I)
                           + np.einsum("jk,il", A2, I)
                          ) 
                - (1.0/35.0)*(
                              np.einsum("ij,kl", I, I)
                            + np.einsum("ik,jl", I, I)
                            + np.einsum("il,jk", I, I)
                            ))
    
    return A4_linear

def quadratic(A2):
    return np.einsum("ij,kl", A2, A2)

    
def hybrid(A2):
    """
    	Literature:
        S. G. Advani and C. L. Tucker,
        “Closure approximations for three‐dimensional structure tensors,” 
        J. Rheol. (N. Y. N. Y)., vol. 34, no. 3, pp. 367–386, 1990, doi: 10.1122/1.550133.
    """
    f = 1.0 - 27.0*np.linalg.det(A2)
    
    A4_hybrid= (1.0 - f)*linear(A2) + f*quadratic(A2)
    
    return A4_hybrid