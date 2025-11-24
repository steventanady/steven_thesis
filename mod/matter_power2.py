# This module is used with camb
from cosmosis.datablock import names, option_section
import numpy as np
from cosmology_mod import Cosmology
import time
import matplotlib.pyplot as plt
from model_independent import model_independent
from model_independent import Growth_Factor

def setup (options) :
    return 0

def linear_evolution(block, a) :
    h0 = block["cosmological_parameters","h0"]
    exponent = block["distances", "exponent"] # The 'n' in the model-indep expf.

    # Read the cheb coeffs from block
    c = []
    i=0
    while True:
        if block.has_value("chebyshev_coefficients", f"c{i}"):
            c.append(block["chebyshev_coefficients", f"c{i}"])
            i += 1
        else:
            break
    c = np.array(c)

    if ("cosmological_parameters", "omh2") in block :
        omega_m = block["cosmological_parameters", "omh2"]/h0**2
    else :
        omega_c = block["cosmological_parameters","omch2"]/h0**2
        omega_b = block["cosmological_parameters","ombh2"]/h0**2
        omega_m = omega_c + omega_b

    # Read the LCDM matter power spectrum from CAMB at the z_max
    pm0 = block["matter_power_lin", "p_k"][-1] 

    # Setting up the tools to evolve matter power model independently
    a_bg = block["distances","a"]
    dplus = Growth_Factor (c,a_bg,h0=h0,om=omega_m,n=exponent) # D+ at each a

    # Calculating model-indep D+
    Dplus = [] # array of D+
    for i in a : 
        Dplus.append(dplus.D_plus(i))
    Dplus = np.array(Dplus)
    Dplus = Dplus/Dplus[0] # Re-Normalize it such that D+(z=6) = 1

    # Evolve matter power model independently
    pm = np.zeros((Dplus.shape[0],pm0.shape[0]))
    for i in range(0, Dplus.shape[0]) :
        pm[i] = pm0*Dplus[i]**2
    return pm

def execute (block, config) :
    # Skip the calculation if the sample gives negative expansion function
    # cause we don't like it (also because the project_2d module crashes)
    if not block["distances", "pos_expf"] :
        return 0
    
    a = block["matter_power_lin","a"]
    pm = linear_evolution(block, a)
    # Skip the calculation if the calculated matter power contains 'nan'
    if np.any(np.isnan(pm)) :
        pm = np.ones_like(pm)
        block["distances", "pos_expf"] = False # Using the same flag to as negative expf to skip, but it doesn't really matter
        print("NAN POWER SPECTRUM ,SKIPPING THE FOLLWING MODULES AND ASSIGNING LOW LIKELIHOOD")
        block.put_grid("matter_power_lin", "z", z, "k_h", k_h, "p_k", pm)
        return 0
    
    z = block["matter_power_lin","z"]
    k_h = block["matter_power_lin","k_h"]

    pm = pm[::-1]
    block.put_grid("matter_power_lin", "z", z, "k_h", k_h, "p_k", pm)

    return 0