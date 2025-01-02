#%% Mods
import sympy as sym 
from fiberOrientationModelling_symbolicComputationTools import symm, generateCCode
from itertools import permutations

sym.init_printing(pretty_print=False)

#%% Funcs
def beta_i(coefDict, Num, II, III):
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


def symA4(A4):
    store = sym.ImmutableDenseNDimArray.zeros(3,3,3,3)
    
    for indices in permutations([0,1,2,3], 4):
        store += sym.permutedims(A4, indices)
        
    return sym.Rational(1, 24)*store
        
#%% Closure
def computeIBOFClosure():
    """
      Literature:      
        D. H. Chung and T. H. Kwon,
        “Invariant-based optimal fitting closure approximation for the numerical prediction of flow-induced fiber orientation”
        J. Rheol. (N. Y. N. Y)., vol. 46, no. 1, pp. 169–194, 2002, 
        doi: 10.1122/1.1423312.
    """

    A2 = symm(sym.Matrix(sym.MatrixSymbol('A', 3, 3)))

    # invariants
    II = sym.Rational(1,2)*(sym.trace(A2)**2 - sym.trace(A2*A2))
    III = A2.det() 
    
    # Stores the coefficients for beta_3, beta_4 and beta_6 in a dictionary. Values are from reference.
    b_coefs={'3':[   2.494090816578600E+01,
                    -4.351011531603290E+02,
                     3.723893356638770E+03,
                     7.034436579164760E+03,
                     8.239951873661060E+05,
                    -1.339319298942450E+05,
                     8.806835153279160E+05,
                    -9.916306907419810E+06,
                    -1.593923962373070E+04,
                     8.009700268497960E+06,
                    -2.370104586892520E+06,
                     3.790105993552670E+07,
                    -3.370108202738210E+07,
                     3.222194162564170E+04,
                    -2.572588058705670E+08,
                     2.144190903444740E+06,
                    -4.492755918514900E+07,
                    -2.131339202233550E+07,
                     1.570767023722040E+09,
                    -2.321534885252980E+04,
                    -3.957693983044730E+09],
              
              '4':[ -4.972177901107540E-01,
                     2.349807975114050E+01,
                    -3.910442513978380E+02,
                     1.539658205935060E+02,
                     1.527729507438190E+05,
                    -2.137552487856460E+03,
                    -4.001389470928120E+03,
                    -1.859493059223080E+06,
                     2.960048652758140E+03,
                     2.477178100543660E+06,
                     1.010139833390620E+05,
                     7.323414942135780E+06,
                    -1.479190276442020E+07,
                    -1.040920721897670E+04,
                    -6.351499296243360E+07,
                    -2.474351062102370E+05,
                    -9.029803789292720E+06,
                     7.249697968073990E+06,
                     4.870934528925950E+08,
                     1.380886909649460E+04,
                    -1.601621786142340E+09],
              
              '6':[  2.341462915709990E+01,
                    -4.120480433725340E+02,
                     3.195532003920890E+03,
                     5.732595943310150E+03,
                    -4.852128030648130E+04,
                    -6.050061135155920E+04,
                    -4.771737400175670E+04,
                     5.990664866898360E+06,
                    -1.106569351765690E+04,
                    -4.605435806806960E+07,
                     2.030429603228740E+06,
                    -5.566061567348350E+07,
                     5.674249110078370E+08,
                     1.289670586862040E+04,
                    -1.527528549565140E+09,
                    -4.993217460925340E+06,
                     1.321248281433330E+08,
                    -1.623599946209830E+09,
                     7.925268498822180E+09,
                     4.667675812929850E+03,
                    -1.280507782794590E+10]
              }
    
    # Calculates beta3, beta4 and beta6
    b_3 = beta_i(b_coefs, '3', II, III)
    b_4 = beta_i(b_coefs, '4', II, III)
    b_6 = beta_i(b_coefs, '6', II, III)
    
    # Calculates beta_1, beta_2 and beta_5 from the 2nd and 3rd invariants of A_ij and from beta_3, beta_4 and beta_6
    b_1 = sym.Rational(3, 5)*(
                                - sym.Rational(1,7)
                                + sym.Rational(1, 5)*b_3*(
                                                              sym.Rational(1, 7) 
                                                            + sym.Rational(4, 7)*II 
                                                            + sym.Rational(8, 3)*III
                                                          )
                                
                                - b_4*( 
                                          sym.Rational(1,  5) 
                                        - sym.Rational(8,  15)*II 
                                        - sym.Rational(14, 15)*III
                                       )
                                
                                - b_6*(
                                          sym.Rational(1,  35) 
                                        - sym.Rational(24, 105)*III 
                                        - sym.Rational(4,  35)*II
                                        + sym.Rational(16, 15)*II*III
                                        + sym.Rational(8,  35)*II**2
                                       )        
                            )
    
    b_2 = sym.Rational(6, 7)*(
                                  1
                                - sym.Rational(1, 5)*b_3*(1 + 4*II)
                                + sym.Rational(7, 5)*b_4*( sym.Rational(1, 6) - II)
                                - b_6*(
                                        - sym.Rational(1, 5)
                                        + sym.Rational(2, 3)*III
                                        + sym.Rational(4, 5)*II
                                        - sym.Rational(8, 5)*II**2
                                      )
                             )
    
    b_5 = (
            - sym.Rational(4, 5)*b_3
            - sym.Rational(7, 5)*b_4
            - sym.Rational(6, 5)*b_6*( 1 - sym.Rational(4, 3)*II)
          )
    
    # 2nd order identity tensor
    I = sym.eye(3)
    
    # 4th order orientation tensor 
    A4_IBOF = (
                   b_1*symA4( sym.tensorproduct(I, I) )
                 + b_2*symA4( sym.tensorproduct(I, A2) )
                 + b_3*symA4( sym.tensorproduct(A2, A2) )
                 + b_4*symA4( sym.tensorproduct(I, (A2*A2) ) )
                 + b_5*symA4( sym.tensorproduct(A2, (A2*A2) ) )
                 + b_6*symA4( sym.tensorproduct((A2*A2), (A2*A2) ) )
              ) 

    return A4_IBOF

#%% Main
def main():

    D = symm(sym.Matrix(sym.MatrixSymbol('D',3,3)))
    
    A4 = computeIBOFClosure()
    
    # Performs the contraction D_kl : A_ijkl 
    D_doubleDot_A4 = sym.tensorcontraction(sym.tensorproduct(D, A4),(0,4),(1,5))
    
    generateCCode(sym.Matrix(D_doubleDot_A4), getOnlyUniqueValues=True)

    return None

#%% Run
if __name__ == "__main__":
    main()
