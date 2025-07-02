from cosmosis.datablock import names, option_section
import numpy as np
import matplotlib.pyplot as plt

def setup(options) :
    return {}

def execute (block, config) :
    block["distances", "pos_expf"] = True
    return 0