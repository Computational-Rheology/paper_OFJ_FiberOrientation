#%% Modules
from .hybrid import hybrid
from .IBOF import IBOF
from .ORE import ORE

#%% Funcs
def create_closure(name):
    
    closures = {"hybrid":   hybrid,
                "IBOF":     IBOF,
                "ORE":      ORE}
    
    if name in closures.keys():
        return closures[name]
    else:
        f_list = "[\n  " + ",\n  ".join(list(closures.keys())) + "\n]"
        raise RuntimeError(f"{name} is not available.\nThe available closure models are:\n{f_list}")
        