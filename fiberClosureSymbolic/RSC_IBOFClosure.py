#%% Mods
import sympy as sym 
from fiberOrientationModelling_symbolicComputationTools import symm, generateCCode
from IBOFClosure import computeIBOFClosure
from rsc import RSC 

#%% Main
def main():   
    
    D = symm(sym.Matrix(sym.MatrixSymbol('D',3,3)))
    
    A4 = RSC(computeIBOFClosure)
    
    D_doubleDot_A4 = sym.tensorcontraction(sym.tensorproduct(D, A4),(0,4),(1,5))
    
    generateCCode(sym.Matrix(D_doubleDot_A4), getOnlyUniqueValues=True)

    return None

#%% Run
if __name__ == "__main__":
    main()

