from cosmosis.datablock import names, option_section
import numpy as np
import matplotlib.pyplot as plt

def setup(options) :
    return {}

def execute (block, config) :
    # flag = block["matter_power_nl","k015"]
    # if flag :
    k = block["matter_power_nl","k_h"]
    pm_nl = block["matter_power_nl","p_k"]
    z = block["matter_power_nl","z"]
    np.savetxt("pk_mycamb_nl_500_zmin05.txt",pm_nl)

    index = np.where(k>0.1)[0]
    pm_nl[:,index]=0
    block.put_grid("matter_power_nl", "z", z, "k_h", k, "p_k", pm_nl)

    # np.savetxt("pk_mycamb_nl_015_zmin05.txt",pm_nl)
    # np.savetxt("k_mycamb_500_zmin05.txt",k)


    return 0