#%% Modules
import numpy as np
from .utils import coefs_fn, fill_symmetry

#%% Funcs
def ORE(A2):
    """
    	Literature:
        B. E. VerWeyst. Numerical Predictions of Flow-Induced Fiber Orientation in 3-D
        Geometries. PhD thesis, University of Illinois at Urbana-Champaign, Urbana, 1998.
    """

    c_coefs={'1':[ 0.636256796880687 ,
                  -1.872662963738140 ,
                  -4.479708731937380 ,
                   11.958956233232000,
                   3.844596924200860 ,
                   11.342092427815900,
                  -10.958262606969100,
                  -20.727799468413200,
                  -2.116232144710040 ,
                  -12.387563285561900,
                   9.815983897167480 ,
                   3.479015105674390 ,
                   11.749291117702600,
                   0.508041387366637 ,
                   4.883665977714890 ],
              
              '2':[ 0.636256796880687 ,
                   -3.315272297421460 ,
                   -3.037099398254060 ,
                    11.827328596885200,
                    6.881539520580440 ,
                    8.436777467783250 ,
                   -15.912066715764100,
                   -15.151587260630700,
                   -6.487289336419260 ,
                   -8.638914192840160 ,
                    9.325203434526610 ,
                    7.746837517132950 ,
                    7.481468706244410 ,
                    2.284765316379580 ,
                    3.597722511342540 ],
              
              '3':[ 2.740532895602530 ,
                   -9.121965097826920 ,
                   -12.257058703625400,
                    34.319901891698700,
                    13.829469912194000,
                    25.868475525388400,
                   -37.702911802938400,
                   -50.275643192748500,
                   -10.880176113317400,
                   -26.963691523971600,
                    27.334679805448800,
                    15.265068614865100,
                    26.113491400537500,
                    3.432138403347790 ,
                    10.611741806606000]
              }
    
    lambda_, e_ = np.linalg.eig(A2)
    
    idx = lambda_.argsort()[::-1]   
    lambda_ = lambda_[idx]
    e_ = e_[:,idx]
            
    store = np.zeros(6)
    store[0] = coefs_fn(c_coefs,'1', lambda_[0], lambda_[1])
    store[1] = coefs_fn(c_coefs,'2', lambda_[0], lambda_[1])
    store[2] = coefs_fn(c_coefs,'3', lambda_[0], lambda_[1])
    
    store[3] = 0.5 * (1 - 2*lambda_[0] + store[0] - store[1] - store[2])
    store[4] = 0.5 * (1 - 2*lambda_[1] - store[0] + store[1] - store[2])
    store[5] = 0.5 * (-1 + 2*lambda_[0] + 2*lambda_[1] - store[0] - store[1] + store[2])
    
    
    
    A4_hat = np.zeros((3,3,3,3))
    for i in range(3):
        A4_hat[i, i, i, i] = store[i]
    
    A_12 =  store[5]
    A_13 =  store[4]
    A_23 =  store[3]

    fill_symmetry(A4_hat, [0,0,1,1], A_12)
    fill_symmetry(A4_hat, [0,0,2,2], A_13)
    fill_symmetry(A4_hat, [1,1,2,2], A_23)
    fill_symmetry(A4_hat, [1,2,2,1], store[3])
    fill_symmetry(A4_hat, [0,2,2,0], store[4])          
    fill_symmetry(A4_hat, [0,1,1,0], store[5])  
    
    
    A4_ORE = np.einsum("im,jn,ko,lp,mnop->ijkl", e_, e_, e_, e_, A4_hat)
    
    return A4_ORE
    
            