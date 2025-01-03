#%% Mods
import sympy as sym 
from fiberOrientationModelling_symbolicComputationTools import symm, generateCCode
from hybridClosure import computeHybridClosure 
from rsc import RSC 

#%% Main
def main():
    
    d = sym.MatrixSymbol('D',3,3)
    
    D = symm(sym.Matrix(d))
    
    A4 = RSC(computeHybridClosure)
    
    D_doubleDot_A4 = sym.tensorcontraction(sym.tensorproduct(D, A4),(0,4),(1,5))
    
    generateCCode(sym.Matrix(D_doubleDot_A4), getOnlyUniqueValues=True)

    return None

#%% Run
if __name__ == "__main__":
    main()

