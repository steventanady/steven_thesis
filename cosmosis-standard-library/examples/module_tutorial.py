from cosmosis.datablock import names, option_section

def setup (options) :
    # This part takes information from the setup parameters of the module in the .ini file
    # Since I am not going to put in any setup parameters, I leave it blank
    return {}

def execute (block, config) :
    # Some values from previous modules
    z = block["distances","z"]
    d_a = block["distances","D_A"]
    H = block["distances","H"]

    # Some constants
    c = 1.0

    # compute something derived from the inputs
    d_v = ((1.0+z)**2 * c * z * d_a**2 / H )**(1./3.)

    block["distances","d_v"] = d_v
    print("d_v[3] = ",d_v[3])

    return 0