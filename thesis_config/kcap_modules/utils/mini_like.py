import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    like_name = options.get_string(option_section, "like_name")
    input_section_name = options.get_string(option_section, "input_section_name", default="likelihood")
    theory60 = options.get_string(option_section, "theory60") # MODDD
        
    return like_name, input_section_name, theory60

def execute(block, config):
    # like_name, input_section_name = config
    like_name, input_section_name, theory60 = config # MODDD, ORIGINALLY ABOVE

    d = block[input_section_name, "data"]
    # d = np.loadtxt(theory60) # MODDD, ORIGINALLY ABOVE, Assigning data as theory from Linear LCDM theta60 result to test if the MCMC works
    mu = block[input_section_name, "theory"]

    inv_cov = block[input_section_name, "inv_covariance"]
    r = d - mu  
    # np.savetxt("d.txt",d)

    chi2 = float(r @ inv_cov @ r)
    if not block["distances", "pos_expf"] :
        chi2 = 9e+100
        print("Since the expansion function/comoving distance is not positive-definite, we ignore this sample by assigning an extremely low likelihood")
    ln_like = -0.5*chi2

    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = ln_like

    return 0

def clean(config):
    pass
