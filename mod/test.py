import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
# a = np.polynomial.Chebyshev([1,2,3], domain=[-1.,110], window=[-1.,1.])
# print(a([1,5]))

z = np.linspace(0.0, 6.0, 200)
a = 1/(1+z)
a = a[::-1]
c1 = 0.5
c2 = 0.5
c3 = 0.5
# a = np.linspace(1/(1+z[-1]),1/(1+z[0]),200)

cheb_poly = np.polynomial.chebyshev.Chebyshev([0.5,0.5,0.5], domain=[a[0],a[-1]])
n = a.shape[0]
h0 = 0.6898
eofa = h0 * cheb_poly(a)**2

d_eofa = cheb_poly.deriv()
d_eofa = h0 * d_eofa(np.linspace(a[0],a[-1],n))

c = 1

    # Calculate comoving distance w(a), which in turn gives a(w) by equalizing the n-th a entry and the n-th w entry.
w = []
for i in range(0,a.shape[0]) :
    x = a[i:]
    w.append(c/h0 * trapezoid(1/(x**2 * eofa[i:]),x))
w = np.array(w)
print(np.all(np.diff(w)<0))
plt.plot(a,w)
plt.savefig("w(a).png")

# print([d_eofa[np.argmin(eofa)-1],d_eofa[np.argmin(eofa)],d_eofa[np.argmin(eofa)+1],d_eofa[np.argmin(eofa)+2]])


# w = []
# for i in range(0,a.shape[0]) :
#     x = a[i:-1]
#     w.append(1/h0 * trapezoid(x,1/(x**2 * eofa[i])))
# w = np.array(w)

# for i in range(0,w.shape[0]) :
#     if(w[i]>1000) :
#         print (a[i]," ",eofa[i]," ",w[i])
# plt.plot(a,eofa)
# plt.ylabel("eofa")
# plt.xlabel("a")
# plt.savefig("test.png")