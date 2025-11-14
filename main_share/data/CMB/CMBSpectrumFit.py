#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=14)  # fontsize of the figure title

plt.style.use('tableau-colorblind10') #clorblind palette for plotting

# import CMB power spectrum from Planck 2018
l_bf, TT_bf = np.loadtxt('COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt', usecols=(0,1), unpack=True)
l, TT, TTm, TTp = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt',unpack=True)

TT_err = (TTm + TTp)/2

# limit spectrum to first peak for fit
l_ = l[(l>50) & (l<350)]
TT_ = TT[(l>50) & (l<350)]
TT_err_ = TT_err[(l>50) & (l<350)]


# define gaussian function
def gaussian(x, a, mu, sigma, offset):
	return a * np.exp( -(x-mu)**2 / (2 * sigma**2)) + offset

# fit gaussian to first peak
guess_gauss = [5700., 220., 100., 600.]
popt, pcov = curve_fit(gaussian, l_, TT_, p0 = guess_gauss, sigma=TT_err_, absolute_sigma=True)


# visualize fit results
# plt.plot(l, TT,'.', label='Planck Data')
plt.errorbar(l,TT,xerr=None, yerr=TT_err, fmt='.', label='Planck Data',
	capsize=0,elinewidth=0.4,markersize=3.0,color='tab:grey',alpha=0.2)
plt.errorbar(l_,TT_,xerr=None, yerr=TT_err_, fmt='.', label='Datasubset for Fit',
	capsize=0,elinewidth=0.4,markersize=3.0,color='tab:grey',alpha=0.5)
plt.plot(l_bf,TT_bf, label='Planck $\Lambda$CDM Fit', color='tab:orange')
plt.plot(l, gaussian(l, *popt), label='Gaussian Fit', color='tab:blue')
plt.grid(True, ls=':', c='#bfbfbf')
plt.xlim(2, 1500)
# plt.xscale('log')
plt.xlabel('Multipole Moment $l$')
plt.ylabel('$\mathcal{D}^{TT}_l$ [$\mu\mathrm{K}^2$]')
plt.title('Power Spectrum of CMB Temperature Map')
plt.legend(loc='upper right').get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.savefig('CMB_fit.pdf', format='pdf')
plt.show()
plt.close()

# calculate distance to CMB in dependance of sound horizon
# d = D_ang/r_d

l0 = popt[1]
l0_err = np.sqrt(pcov[1][1])
sigma = popt[2]

theta = np.pi/l0
theta_err = np.pi*sigma/l0**2

print(l0, l0_err)
results = open("gauss_fit.txt", "w")
results.write('a = '+str(popt[0])+' +- '+str(np.sqrt(pcov[0][0]))+"\n")
results.write('mu = '+str(popt[1])+' +- '+str(np.sqrt(pcov[1][1]))+"\n")
results.write('sigma = '+str(popt[2])+' +- '+str(np.sqrt(pcov[2][2]))+"\n")
results.write('offset = '+str(popt[3])+' +- '+str(np.sqrt(pcov[3][3]))+"\n")
results.write('theta = '+str(theta)+' +- '+str(theta_err)+"\n")

r_d=135.9*1e6
r_d_err=4.8*1e6

z_cmb=1090
a_cmb=1/(1+z_cmb)

mu=5*np.log10(r_d/theta/a_cmb/10)
mu_err=5/np.log(10)*np.sqrt((r_d_err/r_d)**2+(theta_err/theta)**2)

print(theta, theta_err)
print(mu, mu_err)

file = open("cmbSample.d", "w")
file.write(str(z_cmb)+"\t"+str(mu)+"\t"+str(mu_err/2)+"\t"+str(mu_err/2)+"\n")

#Planck 2018 result for theta and corresponding mu

theta = 1.0411/100
theta_err = 0.0003/100

mu=5*np.log10(r_d/theta/a_cmb/10)
mu_err=5/np.log(10)*np.sqrt((r_d_err/r_d)**2+(theta_err/theta)**2)

print(theta, theta_err)
print(mu, mu_err)

file = open("cmbSamplePlanck.d", "w")
file.write(str(z_cmb)+"\t"+str(mu)+"\t"+str(mu_err/2)+"\t"+str(mu_err/2)+"\n")


#Volanthen et al 2010

theta = 0.593*np.pi/180
theta_err = 0.001*np.pi/180

mu=5*np.log10(r_d/theta/a_cmb/10)
mu_err=5/np.log(10)*np.sqrt((r_d_err/r_d)**2+(theta_err/theta)**2)

print(theta, theta_err)
print(mu, mu_err)

file = open("cmbSampleVolanthen.d", "w")
file.write(str(z_cmb)+"\t"+str(mu)+"\t"+str(mu_err/2)+"\t"+str(mu_err/2)+"\n")
