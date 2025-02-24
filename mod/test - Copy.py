import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
import camb
# a = np.polynomial.Chebyshev([1,2,3], domain=[-1.,110], window=[-1.,1.])
# print(a([1,5]))

z = np.linspace(0.0, 6.0, 200)
a = 1/(1+z)
a = a[::-1]
k_h = np.loadtxt("k_h.txt",skiprows=1)
p_k = np.loadtxt("p_k.txt",skiprows=1)
pm = np.loadtxt("pm.txt")
dpl = np.loadtxt("dpl2.txt")

print(pm.shape)
plt.loglog(k_h,p_k[0])
plt.ylabel("P(k,z=0)")
plt.savefig("test2.png")

