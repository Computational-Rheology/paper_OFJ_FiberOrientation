#%% Modules
import numpy as np
import numba as nb
# from itertools import permutations

#%% Funcs
def fill_symmetry(A4, indices, value):
    i,j,k,l = indices
    A4[i,j,k,l] = A4[k,l,i,j] = A4[j,i,k,l] = A4[i,j,l,k] = value
    
def coefs_fn(coefDict, Num, lambda1, lambda2):
    result =  (
                   coefDict[Num][0] 
                 + coefDict[Num][1]*lambda1
                 + coefDict[Num][2]*lambda2
                 + coefDict[Num][3]*lambda1*lambda2
                 + coefDict[Num][4]*(lambda1**2)
                 + coefDict[Num][5]*(lambda2**2)
                 + coefDict[Num][6]*(lambda1**2)*lambda2
                 + coefDict[Num][7]*lambda1*(lambda2**2)
                 + coefDict[Num][8]*(lambda1**3)
                 + coefDict[Num][9]*(lambda2**3)
                 + coefDict[Num][10]*(lambda1**2)*(lambda2**2)
                 + coefDict[Num][11]*(lambda1**3)*lambda2
                 + coefDict[Num][12]*lambda1*(lambda2**3)
                 + coefDict[Num][13]*(lambda1**4)
                 + coefDict[Num][14]*(lambda2**4)
             )
             
    return result

def beta_fn(coefDict, Num, II, III):
    beta =  (
               coefDict[Num][0] 
             + coefDict[Num][1]*II 
             + coefDict[Num][2]*II**2 
             + coefDict[Num][3]*III
             + coefDict[Num][4]*III**2 
             + coefDict[Num][5]*II*III 
             + coefDict[Num][6]*II**2*III 
             + coefDict[Num][7]*II*III**2
             + coefDict[Num][8]*II**3 
             + coefDict[Num][9]*III**3 
             + coefDict[Num][10]*II**3*III 
             + coefDict[Num][11]*II**2*III**2
             + coefDict[Num][12]*II*III**3 
             + coefDict[Num][13]*II**4 
             + coefDict[Num][14]*III**4 
             + coefDict[Num][15]*II**4*III 
             + coefDict[Num][16]*II**3*III**2 
             + coefDict[Num][17]*II**2*III**3 
             + coefDict[Num][18]*II*III**4 
             + coefDict[Num][19]*II**5
             + coefDict[Num][20]*III**5 
          )
    
    return beta

@nb.njit
def symm_A4(A):
    res = np.zeros((3,3,3,3))
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    res[i,j,k,l] = (
                                     A[i, j, k, l] +
                                     A[i, j, l, k] +
                                     A[i, k, j, l] +
                                     A[i, k, l, j] +
                                     A[i, l, j, k] +
                                     A[i, l, k, j] +
                                     A[j, i, k, l] +
                                     A[j, i, l, k] +
                                     A[j, k, i, l] +
                                     A[j, k, l, i] +
                                     A[j, l, i, k] +
                                     A[j, l, k, i] +
                                     A[k, i, j, l] +
                                     A[k, i, l, j] +
                                     A[k, j, i, l] +
                                     A[k, j, l, i] +
                                     A[k, l, i, j] +
                                     A[k, l, j, i] +
                                     A[l, i, j, k] +
                                     A[l, i, k, j] +
                                     A[l, j, i, k] +
                                     A[l, j, k, i] +
                                     A[l, k, i, j] +
                                     A[l, k, j, i]
                                   )
                    
    return (1.0/24)*res       
                    
# list(permutations(["i", "j", "k", "l"]))
