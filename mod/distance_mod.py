from cosmosis.datablock import names, option_section
import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from model_independent import model_independent

def setup (options) :
    # This part takes information from the setup parameters of the module in the .ini file
    zmax = options.get_double(option_section, "zmax")
    # zmid = options.get_double(option_section, "zmid")
    zmin = options.get_double(option_section, "zmin", default = 0.0)    
    nz = options.get_int(option_section, "nz", default = 150)
    nz_mid = options.get_int(option_section, "nz_mid") 
    method = options.get_string(option_section, "method") # Different methods for the model-indep expf. (currently on the inv_cheb works)
    exponent = int(options.get(option_section, "exponent"))
    z = np.linspace(zmin,zmax,nz)
    # z = np.concatenate((np.linspace(zmin, zmid, nz_mid, endpoint = False),
    #                     np.linspace(zmid, zmax, nz-nz_mid)))
    # z = np.loadtxt(options[option_section,"z"])
    a = 1/(1+z)
    a = a[::-1]

    # For the background
    z_bg_max = options.get_double(option_section, "background_zmax")
    z_bg_min = options.get_double(option_section, "background_zmin")  
    nz_bg = options.get_int(option_section, "background_nz")

    z_bg = np.linspace(z_bg_min,z_bg_max,nz_bg)
    a_bg = 1/(1+z_bg) 
    a_bg = a_bg[::-1]

    # input_section_name = options.get_string(option_section, "input_section_name", default="likelihood")
    return {"z":z, "a":a, "z_bg":z_bg, "a_bg":a_bg, "method":method, "exponent":exponent} #, "input_section_name":input_section_name

def execute (block, config) :    
    z = config["z"]
    a = config["a"]
    z_bg = config["z_bg"]
    a_bg = config["a_bg"]
    method = config["method"]
    exponent = config["exponent"]
    h0 = block["cosmological_parameters","h0"]
    if ("cosmological_parameters", "omh2") in block :
        omega_m = block["cosmological_parameters", "omh2"]/h0**2
    else :
        omega_c = block["cosmological_parameters","omch2"]/h0**2
        omega_b = block["cosmological_parameters","ombh2"]/h0**2
        omega_m = omega_c + omega_b
    
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
    
    model = model_independent(c,a_bg,h0,om=omega_m,n=exponent)
    if method == "lcdm" :
        expf = model.lcdm_expf()
    elif method == "pca" :
        expf = model.model_pca(c)
    elif method == "inv_cheb" :
        expf = model.model_expf(c=c)
    
    # Flag for negative expansion function to skip calculation in the following modules
    if np.any(expf<=0) :
        block["distances", "pos_expf"] = False
        print("!!!!!EXPANSION FUNCTION IS NOT POSTIVE-DEFINITE!!!!!")
        print("Assigning a positive-definitive expansion function to let the pipeline run")
        expf = np.linspace(1,2,a_bg.shape[0])
    else :
        block["distances", "pos_expf"] = True
        print("Expansion function is positive-definite, continue as normal")

    c = 299792.458 # km/s

    # Calculate comoving distance w(a), which in turn gives a(w) by equalizing the n-th a entry and the n-th w entry.
    w = model.fast_comoving_distance()
    w = w[::-1]

    # Slow comoving distance
    # w = []
    # for i in range(0,a_bg.shape[0]) :
    #     x = a_bg[i:] #x = a_bg[i:]
    #     w.append(c/(h0*100) * trapezoid(1/(x**2 * expf[i:]),x))
    # w = np.array(w[::-1])
    # print(w)

    # Saving parameters to block
    block["distances", "me"] = True
    block["distances", "d_m"] = w
    block["distances", "a"] = a_bg
    block["distances", "z"] = z_bg
    block["distances", "exponent"] = exponent

    block["matter_power_lin", "a"] = a
    block["matter_power_lin", "z"] = z

    block["cosmological_parameters", "omega_m"] = omega_m

    return 0