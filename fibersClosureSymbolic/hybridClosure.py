#%% Mods
import sympy as sym

from fiberOrientationModelling_symbolicComputationTools import  \
symm, fourthOrderIndexPermutation, generateCCode

sym.init_printing(pretty_print=False)

#%% closure
def computeHybridClosure():
    """      
        Literature:
            S. G. Advani and C. L. Tucker,
            “Closure approximations for three‐dimensional structure tensors,”
            J. Rheol. (N. Y. N. Y)., vol. 34, no. 3, pp. 367–386, 1990, 
            doi: 10.1122/1.550133.

    """

    a = sym.MatrixSymbol('A',3,3)

    A2 = symm(sym.Matrix(a))
        
    I = sym.eye(3)
    
    # Linear Closure
    A4_linear = (
                    sym.Rational(-1,35)*
                    (
                          fourthOrderIndexPermutation('ijkl', sym.tensorproduct(I,I)) 
                        + fourthOrderIndexPermutation('ikjl', sym.tensorproduct(I,I)) 
                        + fourthOrderIndexPermutation('iljk', sym.tensorproduct(I,I))
                    )
                    
                    + sym.Rational(1,7)*
                    (
                          fourthOrderIndexPermutation('ijkl', sym.tensorproduct(A2,I))
                        + fourthOrderIndexPermutation('ikjl', sym.tensorproduct(A2,I))
                        + fourthOrderIndexPermutation('iljk', sym.tensorproduct(A2,I))
                        + fourthOrderIndexPermutation('klij', sym.tensorproduct(A2,I))
                        + fourthOrderIndexPermutation('jlik', sym.tensorproduct(A2,I))
                        + fourthOrderIndexPermutation('jkil', sym.tensorproduct(A2,I))
                    )
                 )
    
    # Quadratic closure
    A4_quadratic =  sym.tensorproduct(A2, A2)
    
    f = 1.0 - 27*A2.det()
    
    A4_hybrid = f*A4_quadratic + (1.0 - f)*A4_linear 
                
    return A4_hybrid
        
#%%  Main
def main():

    d = sym.MatrixSymbol('D',3,3)
    
    D = symm(sym.Matrix(d))
    
    A4 = computeHybridClosure()
    
    # Performs the contraction D_kl : A_ijkl 
    D_doubleDot_A4 = sym.tensorcontraction(sym.tensorproduct(D, A4),(0,4),(1,5))
    
    generateCCode(sym.Matrix(D_doubleDot_A4), getOnlyUniqueValues=True)
      
    return None    
#%% Run        
if __name__ == "__main__":
    main()
