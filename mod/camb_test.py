from cosmosis.datablock import names, option_section
import numpy as np
import matplotlib.pyplot as plt

def setup(options) :
    return

def execute (block, config) :
    p_k = block("matter_power_lin","p_k")
    k = block("matter_power_lin","k_h")
    plt.plot(k,p_k)
    plt.savefig("pofk.png")
