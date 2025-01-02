#%% Mods
import sympy as sym 
import re as re
from sympy.codegen.rewriting import create_expand_pow_optimization

sym.init_printing(pretty_print=False)

#%% Funcs
def symm(T):
    entries = [(0,1), (0,2), (1,2)]

    for entry in entries:
        T[entry[::-1]] = T[entry]
    
    return T


def fourthOrderIndexPermutation(indices, T4):
    res = sym.permutedims(T4, index_order_old=indices, index_order_new='ijkl')
    return res


def place_OF_syntax(lines, vec_symbol, mat_symbol):
    matrix_replaces = { '[0]': '.xx()', '[1]': '.xy()', '[2]': '.xz()',
                        '[3]': '.yx()', '[4]': '.yy()', '[5]': '.yz()',
                        '[6]': '.zx()', '[7]': '.zy()', '[8]': '.zz()' }
                    
    vector_replaces = {'[0]': '.x()', '[1]': '.y()', '[2]': '.z()'}   
    
    for i in range(len(lines)): 
        for vsymb in vec_symbol:
            pattern = fr'{vsymb.name}(\[\d+\])'
            lines[i] = re.sub(pattern,  lambda x: vsymb.name + vector_replaces.get(x.group(1)), lines[i])  
        
        for msymb in mat_symbol:
            pattern = fr'{msymb.name}(\[\d+\])'
            lines[i] = re.sub(pattern,  lambda x: msymb.name + matrix_replaces.get(x.group(1)), lines[i])  
            
    return lines


def generateCCode(Tensor, n_cases=17, OFSyntax=True, getOnlyUniqueValues=True):
        
    sub_exprs, sol = sym.cse(Tensor, sym.numbered_symbols('tmp'))
    if len(sol) != 1:
     raise RuntimeError("something is wrong...")
             
    expand_if_needed = create_expand_pow_optimization(5) 

    lines = ["const scalar " 
             + 
             sym.ccode(expand_if_needed(sym.N(item[1], n=n_cases)), item[0]) 
             for item in sub_exprs]
     
    # Append solution 
    sol_name = "result"
    result = sym.ccode(expand_if_needed(sym.N(sol[0], n=n_cases)), sol_name)
    result = list(filter(None, re.split(";", result)))
    
    # Post processing:        
    free_symbols = Tensor.free_symbols
    mat_symbol = [sym.MatrixSymbol(sol_name, *sol[0].shape)]
    vec_symbol = []
    
    for symb in free_symbols:
        try: 
            if symb.shape == (3,3):
                mat_symbol.append(symb)
            elif symb.shape == (1,3) or symb.shape == (3, 1):
                vec_symbol.append(symb)
        except:
            pass
        
    if OFSyntax:
        lines = place_OF_syntax(lines, vec_symbol, mat_symbol)
        result = place_OF_syntax(result, vec_symbol, mat_symbol)
            
    unique_vals = []
    if getOnlyUniqueValues:
        rhs = []
        for val in result:
            tmp_lhs, tmp_rhs = val.split("=")
            if tmp_rhs not in rhs:
                rhs.append(tmp_rhs)
                unique_vals.append(val + ";")
    else:
        unique_vals = [l + ";" for l in result]
        
    # Print
    for line in lines:
        print(line)
        
    print("\n")
    for val in unique_vals:
        print(val)
    
        
if __name__ == "__main__":
    pass