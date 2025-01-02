#%% Mods
from .folgar_tucker import FT
from .FT_RSC import FT_RSC
from .iARD_RPR import iARD_RPR
from .MRD import MRD
from .MRD_RSC import MRD_RSC
from .pARD_RPR import pARD_RPR


def create_model(name):
    models = {"FT":         FT,
              "FT_RSC":     FT_RSC,
              "iARD_RPR":   iARD_RPR,
              "MRD":        MRD,
              "MRD_RSC":    MRD_RSC,
              "pARD_RPR":   pARD_RPR}
    
    if name in models.keys():
        return models[name]
    else:
        f_list = "[\n  " + ",\n  ".join(list(models.keys())) + "\n]"
        raise RuntimeError(f"{name} is not available.\nThe available models are:\n{f_list}")
        