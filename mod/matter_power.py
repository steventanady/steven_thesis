from cosmosis.datablock import names, option_section
import numpy as np
from scipy.integrate import trapz
from cosmology import Cosmology
import time
import matplotlib.pyplot as plt

def setup (options) :
    # This part takes information from the setup parameters of the module in the .ini file
    
    #############################################
    # Some examples of how this is done:
    # zmax = options[option_section, "zmax"]
	# nz = options[option_section, "nz"]
	# model = options[option_section, "model"]
    #############################################
    # Load matter power spectrum at the furthest source (z = 6)
    pm0 = np.loadtxt(options[option_section,"pm0"])[-1,:]
    k_h = np.loadtxt(options[option_section,"k_h"])
    # k = np.loadtxt("k_h.txt", skiprows=1)

    return {"pm0":pm0, "k_h":k_h}

def linear_evolution(block, config, a) :
    h0 = block["cosmological_parameters","h0"]
    c1 = block["chebyshev_coefficients","c1"]
    c2 = block["chebyshev_coefficients","c2"]
    c3 = block["chebyshev_coefficients","c3"]
    omega_c = block["cosmological_parameters", "omch2"]/h0**2
    omega_b = block["cosmological_parameters", "ombh2"]/h0**2
    omega_m = omega_b + omega_c
    # omega_m = block["cosmological_parameters", "omega_m"]
    pm0 = config["pm0"]
    pm0 = np.array(pm0)
    # print("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
    dplus = Cosmology (c1,c2,c3,a,hubble=h0,omega_matter=omega_m) # D+ at each a

    # Use this for Matthias' code
    # dplus = Cosmology(a[0],hubble = h0, omega_matter = omega_m) 
    Dplus = [] # array of D+
    for i in a : 
        Dplus.append(dplus.D_plus(i))
    Dplus = np.array(Dplus)
    Dplus = Dplus/Dplus[0] # Re-Normalize it such that D+(z=6) = 1
    # np.savetxt("dpl.txt", Dplus)

    pm = np.zeros((Dplus.shape[0],pm0.shape[0]))
    for i in range(0, Dplus.shape[0]) :
        pm[i] = pm0*Dplus[i]**2
    # np.savetxt("pm.txt", pm)
    return pm

# def linear_evolution(growth_factor) :
#     growth_factor = []
#     a = block["distances","a"]
#     h0 = block["cosmological_parameters","h0"]
#     c1 = block["chebyshev_coefficients","c1"]
#     c2 = block["chebyshev_coefficients","c2"]
#     c3 = block["chebyshev_coefficients","c3"]
#     for i in a :
#         Dplus = Cosmology (c1,c2,c3,hubble=h0)
#         growth_factor.append(Dplus.D_plus(i))
#     growth_factor = np.array(growth_factor)

#     # Evolve it
#     pm = config["pm"]
#     for i in a :
#         pm_evolved = pm * growth_factor(i)
    
#     return pm_evolved

def execute (block, config) :
    #############################################
    # Examples :
    # a = config["a"] #config is used to get values from setup

    # Some values from previous modules
    # z = block["distances","z"]
    # d_a = block["distances","D_A"]
    # H = block["distances","H"]

    # Save to data block
    # d_v = ((1.0+z)**2 * c * z * d_a**2 / H )**(1./3.)
    # block["distances","d_v"] = d_v
    #############################################
    a = block["matter_power_lin","a"]
    pm = linear_evolution(block, config, a)
    pm = pm[::-1]
    # np.savetxt("pm_lin_me.txt",pm)

    # block["matter_power_lin","p_k"] = pm
    z = block["matter_power_lin","z"]
    k_h = config["k_h"]

    block.put_grid("matter_power_lin", "z", z, "k_h", k_h, "p_k", pm)

    # Non-linear part, simply for the sake to run the pipeline
    pm_nl = pm
    pm_nl = np.zeros_like(pm)
    block.put_grid("matter_power_nl", "z", z, "k_h", k_h, "p_k", pm_nl)

    return 0