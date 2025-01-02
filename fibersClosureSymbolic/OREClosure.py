#%% Mods
import sympy as sym
from fiberOrientationModelling_symbolicComputationTools import symm, generateCCode

sym.init_printing(pretty_print=False)

def fillInSymmetry(A4, indices, value):
    i,j,k,l = indices
    A4[i,j,k,l] = A4[k,l,i,j] = A4[j,i,k,l] = A4[i,j,l,k] = value
    
    
def Coefs(coefDict, Num, lambda1, lambda2):
    res =  (
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
             
    return res

#%% Closure
def computeOREClosure():
    """
        Literature:
            B. E. VerWeyst. Numerical Predictions of Flow-Induced Fiber Orientation in 3-D
            Geometries. PhD thesis, University of Illinois at Urbana-Champaign, Urbana,
            IL, 1998.
    """
    
    eigVec = sym.Matrix(sym.MatrixSymbol('eigenVector',3,3))
    
    eigVal = sym.Matrix(sym.MatrixSymbol('eigenValue',3,1))

    # Put variables in OpenFoam order
    tmpEigVal = eigVal.copy()
    
    # sort eigenValues
    eigVal[0] = tmpEigVal[2]
    eigVal[1] = tmpEigVal[1]
    eigVal[2] = tmpEigVal[0]
    
    # sort eigenVectors
    eigVec = eigVec.T
    tmpEigVec = eigVec.copy()
    eigVec[:, 0] = tmpEigVec[:, 2]
    eigVec[:, 1] = tmpEigVec[:, 1]
    eigVec[:, 2] = tmpEigVec[:, 0]
    


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
    
    
    
    # Get the 3 independent entries    
    A_11 = Coefs(c_coefs, '1', eigVal[0], eigVal[1])
    A_22 = Coefs(c_coefs, '2', eigVal[0], eigVal[1])
    A_33 = Coefs(c_coefs, '3', eigVal[0], eigVal[1])
    
    A_44 = sym.Rational(1,2)*(1 - 2*eigVal[0] + A_11 - A_22 - A_33)
    A_55 = sym.Rational(1,2)*(1 - 2*eigVal[1] - A_11 + A_22 - A_33)
    A_66 = sym.Rational(1,2)*(-1 + 2*eigVal[0] + 2*eigVal[1] - A_11 - A_22 + A_33)

    A_12 = A_21 = A_66
    A_13 = A_31 = A_55
    A_23 = A_32 = A_44
    
    A4_hat = sym.MutableDenseNDimArray.zeros(3,3,3,3) 
    A4_hat[0,0,0,0] = A_11
    A4_hat[1,1,1,1] = A_22
    A4_hat[2,2,2,2] = A_33
      
    fillInSymmetry(A4_hat, [0,0,1,1], A_12)
    fillInSymmetry(A4_hat, [0,0,2,2], A_13)
    fillInSymmetry(A4_hat, [1,1,2,2], A_23)
    fillInSymmetry(A4_hat, [1,2,2,1], A_44)
    fillInSymmetry(A4_hat, [0,2,2,0], A_55)          
    fillInSymmetry(A4_hat, [0,1,1,0], A_66)   
       
    A4_ORE = sym.MutableDenseNDimArray.zeros(3,3,3,3) 
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for o in range(3): 
                                for p in range(3):
                                    A4_ORE[i,j,k,l] += eigVec[i,m]*eigVec[j,n]*eigVec[k,o]*eigVec[l,p]*A4_hat[m,n,o,p]
                                    
                     
    return A4_ORE.as_immutable()
        
# Main                
def main():
    
    D = symm(sym.Matrix(sym.MatrixSymbol('D',3,3)))
    
    A4 = computeOREClosure()
    
    # Performs the contraction D_kl : A_ijkl 
    D_doubleDot_A4 = sym.tensorcontraction(sym.tensorproduct(D, A4),(0,4),(1,5))
    
    generateCCode(sym.Matrix(D_doubleDot_A4), getOnlyUniqueValues=True)
      
    return None
        
#%% Run        
if __name__ == "__main__":
    main()

