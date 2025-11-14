from cosmosis.datablock import names, option_section
import numpy as np
import matplotlib.pyplot as plt
from cosmology_mod import Cosmology

def setup(options) :
    return {}

def execute (block, config) :

    # Comparing the data and theory result of the quantity used to test likelihood
    # data = block["scale_cuts_output","data"]
    # theory = block["scale_cuts_output","theory"]
    # cov = block["scale_cuts_output", "covariance"]
    # diff = data - theory
    # diff_frac = diff/data
    # diff_frac2 = diff/theory
    # output = np.column_stack((data,theory,diff,diff_frac,diff_frac2))
    # np.savetxt("cov60.txt", cov)
    # np.savetxt("data-theory_SN_BAO_QSO_Marie.txt",output, header = "data theory diff diff_frac diff_frac2")
    # np.savetxt("data_theta60.txt", data)
    # np.savetxt("theory60.txt", theory)

    # Generate output for power spectrum
    pm_lin = block["matter_power_lin","p_k"]
    np.savetxt("pm_lin_kids2.txt",pm_lin)
    # pm_nl = block["matter_power_nl","p_k"]
    # np.savetxt("pm.txt",pm_nl)

    # Save matter power from camb
    # pm_lin = block["matter_power_lin","p_k"]
    # z = block["matter_power_lin","z"]
    # k = block["matter_power_lin","k_h"]
    # np.savetxt("z_iir.txt", z)
    # np.savetxt("p_iir_pre_ext.txt",pm_lin)
    # np.savetxt("k_pre_ext.txt",k)

    # Save output from projection
    # for i in range(1,6) :
    #     for j in range(1,i+1) :
    #         name = f"bin_{i}_{j}"
    #         c_ell = block["shear_cl",name]
    #         np.savetxt(f"c_ell_lcdm_{i}_{j}.txt",c_ell)
    # c_ell = block["shear_cl", "bin_5_5"]
    # np.savetxt("c_ell_kids_lin55.txt",c_ell)

    # TEST
    # omch2 = block["cosmological_parameters","omch2"]
    # ombh2 = block["cosmological_parameters","ombh2"]
    # h0 = block["cosmological_parameters","h0"]
    # print ("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ omch2 + ombh2= ", omch2+ombh2)
    # print ("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ omega_m = ", block["cosmological_parameters","omega_m"]*h0**2)

    # zm = block["matter_power_lin","z"]
    # k = block["matter_power_lin","k_h"]
    # zd = block["distances","z"]
    # a = block["distances","a"]

    # np.savetxt("z_matter015.txt",zm)
    # np.savetxt("z_distance.txt",zd)
    # np.savetxt("a_distance.txt",a)
    # np.savetxt("k_h1.0.txt",k)

    # Save comoving distance
    # dm = block["distances","d_m"]
    # np.savetxt("w_kids.txt",dm)

    # n_cheb = 4
    # c = np.zeros(n_cheb)
    # for i in range(n_cheb) :
    #     c[i] = block["chebyshev_coefficients",f"c{i}"]
    # a = block["distances", "a"]
    # expf = Cosmology(c,a).expansion_function(a)
    # np.savetxt("expf.txt",expf)

    return 0