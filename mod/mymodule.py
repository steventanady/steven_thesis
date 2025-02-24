from cosmosis.datablock import names, option_section
import numpy as np
from scipy.integrate import trapz

def setup (options) :
    # This part takes information from the setup parameters of the module in the .ini file
    # Since I am not going to put in any setup parameters, I leave it blank
    
    #############################################
    # Some examples of how this is done:
    # zmax = options[option_section, "zmax"]
	# nz = options[option_section, "nz"]
	# model = options[option_section, "model"]
    #############################################

    zmax = options[option_section, "zmax"]
    nz = options[option_section, "nz"]
    z = np.linspace(0,zmax,nz)
    a = 1/(1+z)
    return {"zmax":zmax, "nz":nz, "z":z, "a":a}



def execute (block, config) :
    #############################################
    # Examples :
    # Some values from previous modules
    # z = block["distances","z"]
    # d_a = block["distances","D_A"]
    # H = block["distances","H"]

    # Save to data block
    # d_v = ((1.0+z)**2 * c * z * d_a**2 / H )**(1./3.)
    # block["distances","d_v"] = d_v
    #############################################

    a = config["a"] #config is used to get values from setup
    zmax = config["zmax"]
    nz = config["nz"]
    z = config["z"]
    h0 = block["cosmological_parameters","h0"]

    c1 = block["chebyshev_coefficients","c1"]
    c2 = block["chebyshev_coefficients","c2"]
    c3 = block["chebyshev_coefficients","c3"]

    expf_of_a = h0*np.polynomial.Chebyshev([c1,c2,c3], domain:=[0,a[-1]])
    c = 1

    # Calculate w(a), which in turn gives a(w) by equalizing the n-th a entry and the n-th w entry.
    w = []
    for i in range(0,a.shape[0]) :
        x = a[0:i]
        w.append(c/h0 * trapz(x,1/(x**2 * expf_of_a(x))))
    w = np.array(w) 



    print("Great!")
    return 0