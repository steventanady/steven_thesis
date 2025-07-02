import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.optimize import curve_fit
import scipy.interpolate as interp
import camb
import sys
sys.path.append("/home/steven/thesis/mod")
from model_independent import model_independent
from model_independent import Growth_Factor

model = model_independent([1,2,3],np.array(0.15,0.5,1))
test = model.model_expf()