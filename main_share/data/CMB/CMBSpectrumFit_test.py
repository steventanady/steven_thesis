#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# import CMB power spectrum from Planck 2015
l, TT, TE, EE, BB, PP = np.loadtxt('COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt', unpack=True)

# visualize power spectrum

l_ = l[10:405]
TT_ = TT[10:405]

plt.plot(l, TT,'.', markersize=1)
plt.plot(l_, TT_,'.', markersize=1)
plt.xlabel(u'multipole moment $l$', fontsize=9)
plt.ylabel(u'Power Spectrum of Temperature Map $C^{TT}_l$', fontsize=9)
#plt.legend(loc='upper left', fontsize=8)
#plt.savefig('NAME.pdf', format='pdf')
#plt.show()
plt.close()

# define gaussian function

def gaussian5(x, a1, mu1, sigma1, a2, mu2, sigma2, a3, mu3, sigma3, a4, mu4, sigma4, a5, mu5, sigma5, offset):
	return a1 * np.exp( -(x-mu1)**2 / (2 * sigma1**2)) + \
 		a2 * np.exp( -(x-mu2)**2 / (2 * sigma2**2)) + \
		a3 * np.exp( -(x-mu3)**2 / (2 * sigma3**2)) + \
		a4 * np.exp( -(x-mu4)**2 / (2 * sigma4**2)) + \
		a5 * np.exp( -(x-mu5)**2 / (2 * sigma5**2)) + offset

def gaussian(x, a, mu, sigma, offset):
	return a * np.exp( -(x-mu)**2 / (2 * sigma**2)) + offset


# fit gaussian to first peak

guess_gauss5 = [5300., 220., 90., 2200., 540., 100., 2100., 810., 100., 800., 1130., 150., 400., 1400., 150., 50.] 
guess_gauss = [5700., 220., 100., 600.]

popt5, pcov5 = curve_fit(gaussian5, l, TT, p0 = guess_gauss5)
popt, pcov = curve_fit(gaussian, l_, TT_, p0 = guess_gauss)


# visualize fit results

plt.plot(l, TT,'.', markersize=1, label='data')
#plt.plot(l_, TT_,'.', markersize=1, label='data')
plt.plot(l, gaussian5(l, *popt5), label='n=5 gaussian fit', ls=':', c='r')
plt.plot(l, gaussian(l, *popt), label='gaussian fit', ls='--', c='g')
plt.xlabel(u'multipole moment $l$', fontsize=9)
plt.ylabel(u'Power Spectrum of Temperature Map $C^{TT}_l$', fontsize=9)
plt.legend(loc='upper right', fontsize=8)
plt.savefig('CMBSpectrum_GaussianFit.pdf', format='pdf')
plt.show()
plt.close()
