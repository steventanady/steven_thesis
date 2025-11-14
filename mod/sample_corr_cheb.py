from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
    cov_file = options.get_string(option_section,"cov_file")
    cov = np.loadtxt(cov_file)
    return cov

def execute(block,config):
    cov = config

    # Load the chebyshev coefficients
    c=[]
    i=0
    while True:
        if block.has_value("uncorr_chebyshev_coefficients", f"c{i}"):
            c.append(block["uncorr_chebyshev_coefficients", f"c{i}"])
            i += 1
        else:
            break
    c = np.array(c)

    # Perform Cholesky decomposition
    L = np.linalg.cholesky(cov)
    c_corr = L @ c
    
    # Write correlated coefficients to block
    for i in range(c_corr.shape[0]):
        block["chebyshev_coefficients",f"c{i}"] = c_corr[i]

    return 0