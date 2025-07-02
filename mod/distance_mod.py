from cosmosis.datablock import names, option_section
import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from model_independent import model_independent

def setup (options) :
    # This part takes information from the setup parameters of the module in the .ini file
    zmax = options.get_double(option_section, "zmax", default = 3.0)
    zmid = options.get_double(option_section, "zmid")
    zmin = options.get_double(option_section, "zmin", default = 0.0)    
    nz = options.get_int(option_section, "nz", default = 150)
    nz_mid = options.get_int(option_section, "nz_mid") 
    method = options.get_string(option_section, "method")
    # z = np.linspace(zmin,zmax,nz)
    z = np.concatenate((np.linspace(zmin, zmid, nz_mid, endpoint = False),
                        np.linspace(zmid, zmax, nz-nz_mid)))
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
    return {"z":z, "a":a, "z_bg":z_bg, "a_bg":a_bg, "method":method} #, "input_section_name":input_section_name

def execute (block, config) :    
    z = config["z"]
    a = config["a"]
    z_bg = config["z_bg"]
    a_bg = config["a_bg"]
    method = config["method"]
    h0 = block["cosmological_parameters","h0"]
    omega_c = block["cosmological_parameters","omch2"]/h0**2
    omega_b = block["cosmological_parameters","ombh2"]/h0**2
    omega_m = omega_c + omega_b
    
    c = []
    i=0
    while True:
        if block.has_value("chebyshev_coefficients", f"c{i}"):
            c.append(block["chebyshev_coefficients", f"c{i}"])
            i += 1
        else:
            break
    c = np.array(c)
    
    model = model_independent(c,a_bg,h0,om=omega_m)
    if method == "lcdm" :
        expf = model.lcdm_expf()
    elif method == "pca" :
        expf = model.model_pca(c)
    elif method == "inv_cheb" :
        expf = model.model_expf(c=c)
    np.savetxt('expf.txt',expf)
    if np.any(expf<=0) :
        block["distances", "pos_expf"] = False
        print("!!!!!EXPANSION FUNCTION IS NOT POSTIVE-DEFINITE!!!!!")
        print("Assigning a positive-definitive expansion function to let the pipeline run")
        expf = np.linspace(1,2,a_bg.shape[0])
    else :
        block["distances", "pos_expf"] = True
        print("Expansion function is positive-definite, continue as normal")
    

    ##################################################################################################################
    # # LCDM
    # omega_r = 8.5e-5
    # omega_c = block["cosmological_parameters","omch2"]/h0**2
    # omega_b = block["cosmological_parameters","ombh2"]/h0**2
    # omega_m = omega_c + omega_b
    # # omega_m = block["cosmological_parameters", "omega_m"]
    # omega_l = 1.0-omega_m # 0.6889
    # omega_c = 0 # 1.0-omega_r-omega_m-omega_l
    # # expf = np.sqrt (omega_r/a_bg**4 + omega_m/a_bg**3 + omega_l + omega_c/a_bg**2)
    # # np.savetxt("lcdm_expf.txt",expf)

    # ## Model-independent
    # # cheb_poly = np.polynomial.chebyshev.Chebyshev([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15], domain=[a[0],a[-1]])
    # # cheb_poly = np.polynomial.chebyshev.Chebyshev([c1,c2,c3,c4,c5], domain=[0,1]) #domain=[a_bg[0],a_bg[-1]]
    
    # n_cheb = 4
    # c = np.zeros(n_cheb)
    # for i in range(n_cheb) :
    #     c[i] = block["chebyshev_coefficients",f"c{i}"]
    # def cheb(c,a):
    #     x = np.arccos(a)
    #     result = 0
    #     for i in range(len(c)) :
    #         if i ==0 :
    #             result += (1/np.pi)**(1/2)*c[i]
    #         else :
    #             result += (2/np.pi)**(1/2)*c[i]*np.cos(i*x)
    #     return result

    # da = (1-a_bg[0])/a_bg[0]
    # x = (a_bg-a_bg[0])/(a_bg[-1]-a_bg[0])
    # expf = da*h0/(1+da*x)**2/(cheb(c,2*x-1))
    # expf = expf/expf[-1]
    # # np.savetxt("model_expf.txt",expf)
    # # np.savetxt("a_bg.txt",a_bg)

    # #Check for whether the expansion function is positive-definite
    # if np.any(expf<=0) :
    #     block["distances", "pos_expf"] = False
    #     print("!!!!!EXPANSION FUNCTION IS NOT POSTIVE-DEFINITE!!!!!")
    #     print("Assigning a positive-definitive expansion function to let the pipeline run")
    #     expf = np.linspace(1,2,a_bg.shape[0])
    # else :
    #     block["distances", "pos_expf"] = True
    #     print("Expansion function is positive-definite, continue as normal")
    #     # expf = h0 * cheb_poly(a)    
    ##################################################################################################################

    c = 299792.458 # km/s
    # dh = c/100 # 

    # Calculate comoving distance w(a), which in turn gives a(w) by equalizing the n-th a entry and the n-th w entry.
    w = []
    for i in range(0,a_bg.shape[0]) :
        x = a_bg[i:] #x = a_bg[i:]
        w.append(c/(h0*100) * trapz(1/(x**2 * expf[i:]),x))
    w = np.array(w[::-1])

    # for i in range(0,a_bg.shape[0]) :
    #     x = a_bg[0:i] #x = a_bg[i:]
    #     w.append(-c/(h0*100) * trapz(1/(x**2 * expf[0:i]),x))

    block["distances", "me"] = True
    block["distances", "d_m"] = w
    block["distances", "a"] = a_bg
    block["distances", "z"] = z_bg

    block["matter_power_lin", "a"] = a
    block["matter_power_lin", "z"] = z

    block["cosmological_parameters", "omega_m"] = omega_m

    return 0