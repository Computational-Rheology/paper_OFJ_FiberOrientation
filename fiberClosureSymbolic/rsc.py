#%% Mods
import sympy as sym 
sym.init_printing(pretty_print=False)

#%% RSC
def RSC(closureFunc):
    # Creates a symbolic scalar for the slow-down parameter
    k=sym.Symbol('k')
    
    # Creates a Symbolic Matrix to store eigenVectors
    e = sym.MatrixSymbol('eigenVector',3,3)
    eigVec = sym.Matrix(e)
    
    # Creates a Symbolic Matrix to store eigenValues
    w = sym.MatrixSymbol('eigenValue',3,1)
    eigVal = sym.Matrix(w)
    
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
    
    
    A4 = closureFunc()
    
    # Calculate M and L
    M = sym.MutableDenseNDimArray.zeros(3,3,3,3)
    L = sym.MutableDenseNDimArray.zeros(3,3,3,3)
    
    for i in range(3):
        e_vector =  sym.Array(sym.flatten(eigVec[:, i]))
        M += sym.tensorproduct(e_vector, e_vector, e_vector, e_vector)
        L += eigVal[i]*sym.tensorproduct(e_vector, e_vector, e_vector, e_vector)
    
    tmp1 = sym.tensorcontraction(sym.tensorproduct(M, A4), (2,4), (3,5))

    tensorCombination = A4 + (1 - k)*(L - tmp1)

    return tensorCombination

if __name__ == "__main__":
    pass

                







