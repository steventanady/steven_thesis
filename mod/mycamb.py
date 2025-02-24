from cosmosis.datablock import names, option_section
import numpy as np
import matplotlib.pyplot as plt
import camb
from camb import model, initialpower

def setup (options) :
    # This part takes information from the setup parameters of the module in the .ini file
    
    #############################################
    # Some examples of how this is done:
    # zmax = options[option_section, "zmax"]
	# nz = options[option_section, "nz"]
	# model = options[option_section, "model"]
    #############################################
    zmax = options.get_double(option_section, "zmax", default = 3.0)
    zmid = options.get_double(option_section, "zmid")
    zmin = options.get_double(option_section, "zmin", default = 0.0)    
    nz = options.get_int(option_section, "nz", default = 150)
    nz_mid = options.get_int(option_section, "nz_mid")
    z = np.concatenate((np.linspace(zmin, zmid, nz_mid, endpoint = False),
                        np.linspace(zmid, zmax, nz-nz_mid)))
    kmax = options.get_double(option_section, "kmax")
    return {"z":z, "kmax":kmax}

def execute(block, config) :
    h0 = block["cosmological_parameters", "h0"]
    omch2 = block["cosmological_parameters", "omch2"]
    ombh2 = block["cosmological_parameters", "ombh2"]
    ns = block["cosmological_parameters", "n_s"]
    z = config["z"]


    H0 = h0*100

    pars1 =  camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, ns=ns)
    #Note non-linear corrections couples to smaller scales than you want
    kmax = config["kmax"]
    pars1.set_matter_power(redshifts=z, kmax=kmax)

    # #Linear spectra
    # pars1.NonLinear = model.NonLinear_none
    results = camb.get_results(pars1)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=kmax, npoints = 200)
    block.put_grid("matter_power_lin", "z", z, "k_h", kh, "p_k", pk)
    # s8 = np.array(results.get_sigma8())

    #Non-Linear spectra (Halofit)
    pars1.NonLinear = model.NonLinear_both
    results.calc_power_spectra(pars1)
    kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=kmax, npoints = 200)

    block.put_grid("matter_power_nl", "z", z_nonlin, "k_h", kh_nonlin, "p_k", pk_nonlin)

    return 0