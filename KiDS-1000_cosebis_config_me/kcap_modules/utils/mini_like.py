import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    like_name = options.get_string(option_section, "like_name")
    input_section_name = options.get_string(option_section, "input_section_name", default="likelihood")
    # theory60 = options.get_string(option_section, "theory60") # MODDD
        
    return like_name, input_section_name

def execute(block, config):
    like_name, input_section_name = config
    # like_name, input_section_name, theory60 = config # MODDD, ORIGINALLY ABOVE

    # Skip if expansion function is not positive-definite
    # if not block("distances", "pos_expf") :
    #     print("Assigning extremely small likelihood")
    #     block[names.data_vector, like_name+"_CHI2"] = 
    #     block[names.likelihoods, like_name+"_LIKE"] = 
    #     return 0
    d = block[input_section_name, "data"]
    # d = np.loadtxt(theory60) # MODDD, ORIGINALLY ABOVE
    mu = block[input_section_name, "theory"]

    inv_cov = block[input_section_name, "inv_covariance"]
    r = d - mu

    chi2 = float(r @ inv_cov @ r)
    # chi2 = float(r @ r) # MODDD
    ln_like = -0.5*chi2

    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = ln_like

    return 0

def clean(config):
    pass
