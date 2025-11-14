from cosmosis.datablock import names, option_section
import numpy as np
from cosmology_mod import Cosmology
import time
import matplotlib.pyplot as plt
from model_independent import model_independent
from model_independent import Growth_Factor
def setup (options) :
    pm0 = np.loadtxt(options[option_section,"pm0"])
    k_h = np.loadtxt(options[option_section,"k_h"])

    return {"pm0":pm0, "k_h":k_h}

def linear_evolution(block, config, a) :
    h0 = block["cosmological_parameters","h0"]
    exponent = block["distances", "exponent"]
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

    # Use the following when I am applying the model-independent analysis
    pm0 = config["pm0"][-1,:]
    pm0 = np.array(pm0)

    # zz = np.linspace(0,1100,25000)
    # aa = 1/(1+zz)
    # aa = aa[::-1]
    dplus = Growth_Factor (c,a,h0=h0,om=omega_m,n=exponent)
    
    Dplus = [] # array of D+
    for i in a : 
        Dplus.append(dplus.D_plus(i))
    Dplus = np.array(Dplus)
    Dplus = Dplus/Dplus[0] # Re-Normalize it such that D+(z=6) = 1
    # np.savetxt("dpl_many.txt", Dplus)

    pm = np.zeros((Dplus.shape[0],pm0.shape[0]))
    for i in range(0, Dplus.shape[0]) :
        pm[i] = pm0*Dplus[i]**2
    np.savetxt("pm.txt", pm)
    print(Dplus)

    return pm

def execute (block, config) :
    # Skip if expansion function is not positive-definite
    
    if not block["distances", "pos_expf"] :
        pm0 = config["pm0"]
        pm0 = np.array(pm0)
        pm = np.zeros_like(pm0)
        z = block["matter_power_lin","z"]
        k_h = config["k_h"]
        block.put_grid("matter_power_lin", "z", z, "k_h", k_h, "p_k", pm)
        return 0

    a = block["matter_power_lin","a"]
    pm = linear_evolution(block, config, a)

    # block["matter_power_lin","p_k"] = pm
    z = block["matter_power_lin","z"]
    k_h = config["k_h"]

    # Skip if matter power spectrum gives nan
    if np.any(np.isnan(pm)) :
        pm = np.ones_like(pm)
        block["distances", "pos_expf"] = False
        print("NAN POWER SPECTRUM")
        block.put_grid("matter_power_lin", "z", z, "k_h", k_h, "p_k", pm)
        return 0
    
    pm = pm[::-1]
    # np.savetxt("pkz2.txt",pm)
    block.put_grid("matter_power_lin", "z", z, "k_h", k_h, "p_k", pm)

    # Non-linear part, simply for the sake to run the pipeline
    # pm_nl = pm
    # pm_nl = np.zeros_like(pm)
    # block.put_grid("matter_power_nl", "z", z, "k_h", k_h, "p_k", pm_nl)

    return 0