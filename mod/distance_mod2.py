# This module is used with camb
from cosmosis.datablock import names, option_section
import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt

def setup (options) :
    # This part takes information from the setup parameters of the module in the .ini file
    
    #############################################
    # Some examples of how this is done:
    # zmax = options[option_section, "zmax"]
	# nz = options[option_section, "nz"]
	# model = options[option_section, "model"]
    #############################################

    # zmax = options.get_double(option_section, "zmax", default = 3.0)
    # zmid = options.get_double(option_section, "zmid")
    # zmin = options.get_double(option_section, "zmin", default = 0.0)    
    # nz = options.get_int(option_section, "nz", default = 150)
    # nz_mid = options.get_int(option_section, "nz_mid") 
    # z = np.linspace(zmin,zmax,nz)
    # z = np.concatenate((np.linspace(zmin, zmid, nz_mid, endpoint = False),
                        # np.linspace(zmid, zmax, nz-nz_mid)))
    # z = np.loadtxt(options[option_section,"z"])
    # a = 1/(1+z)
    # a = a[::-1]

    # For the background
    # z_bg_max = options.get_double(option_section, "background_zmax")
    # z_bg_min = options.get_double(option_section, "background_zmin")  
    # nz_bg = options.get_int(option_section, "background_nz")

    # z_bg = np.linspace(z_bg_min,z_bg_max,nz_bg)
    # a_bg = 1/(1+z_bg) 
    # a_bg = a_bg[::-1]

    # input_section_name = options.get_string(option_section, "input_section_name", default="likelihood")
    return {} #, "input_section_name":input_section_name

def execute (block, config) :

    # d = block["likelihood", "data"]
    # np.savetxt("dataaaaaa.txt",data)

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

    # z = config["z"]
    # a = config["a"]
    # z_bg = config["z_bg"]
    # a_bg = config["a_bg"]
    z_bg = block["distances","z"]
    a_bg = 1/(1+z_bg)
    a_bg = a_bg[::-1]

    h0 = block["cosmological_parameters","h0"]
    c1 = block["chebyshev_coefficients","c1"]
    c2 = block["chebyshev_coefficients","c2"]
    c3 = block["chebyshev_coefficients","c3"]

    # LCDM
    omega_r = 8.5e-5
    omega_c = block["cosmological_parameters","omch2"]/h0**2
    omega_b = block["cosmological_parameters","ombh2"]/h0**2
    omega_m = omega_c + omega_b
    # omega_m = block["cosmological_parameters", "omega_m"]
    omega_l = 1.0-omega_m # 0.6889
    omega_c = 0 # 1.0-omega_r-omega_m-omega_l
    eofa = np.sqrt (omega_r/a_bg**4 + omega_m/a_bg**3 + omega_l + omega_c/a_bg**2)

    

    # Model-independent
    # cheb_poly = np.polynomial.chebyshev.Chebyshev([c1,c2,c3], domain=[a[0],a[-1]])
    # n = a.shape[0]
    # eofa = h0 * cheb_poly(a)**2
    c = 299792.458 # km/s
    # dh = c/100 # 

    # Calculate comoving distance w(a), which in turn gives a(w) by equalizing the n-th a entry and the n-th w entry.
    w = []
    for i in range(0,a_bg.shape[0]) :
        x = a_bg[i:] #x = a_bg[i:]
        w.append(c/(h0*100) * trapz(1/(x**2 * eofa[i:]),x))
    w = np.array(w[::-1])

    # for i in range(0,a_bg.shape[0]) :
    #     x = a_bg[0:i] #x = a_bg[i:]
    #     w.append(-c/(h0*100) * trapz(1/(x**2 * eofa[0:i]),x))

    # block["distances", "me"] = True

    # test = block["distances","D_M"]
    # print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ",test[-1])
    # test = block["distances","d_m"]
    # print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ",test[-1])
    block["distances", "d_m"] = w # d_m or D_M is fine, it's not case-sensitive
    # test = block["distances","d_m"]
    # print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ",test[-1])
    # test = block["distances","D_M"]
    # print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ",test[-1])
    
    # block["distances", "a"] = a_bg
    # block["distances", "z"] = z_bg

    # block["matter_power_lin", "a"] = a
    # block["matter_power_lin", "z"] = z

    block["cosmological_parameters", "omega_m"] = omega_m # This is somehow affecting the final result

    return 0