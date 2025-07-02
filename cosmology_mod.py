## Import necessary packages

from scipy.integrate import quad
from scipy.integrate import odeint
import numpy as np

## Class providing functions essential for standard, Friedmann-Lemaître
#  cosmological models
#
#  Matthias Bartelmann, U. Heidelberg, 2023

class Cosmology:
    ## Constructor with default arguments setting the density parameters and
    #  the Hubble constant (in units of 100 km/s/Mpc to the values published
    #  by the Planck project in 2018
    #
    #  omega_matter is the sum of the baryonic and dark matter densities
    def __init__ (self, \
        omega_radiation = 8.5e-5, \
        omega_matter = 0.3111, \
        omega_lambda = 0.6889, \
        hubble = 0.6766):
        self.omega_r  = omega_radiation
        self.omega_m  = omega_matter
        self.omega_l  = omega_lambda
        self.omega_b  = 0.04
        self.omega_c  = 1.0-self.omega_r-self.omega_m-self.omega_l
        self.hubble   = hubble
        self.rho_mean = 2.7594e11*self.hubble**2*self.omega_m
        self.a_equal  = self.omega_r/self.omega_m
        self.n_growth = self.growth_index ()
        self.a_ini    = 1.0e-3
        self.normalize ()
    
    ## Normalization routine setting initial values for the growth factor,
    #  the effective particle mass, and the expansion function
    def normalize (self):
        self.D_plus_0 = 1.0
        self.D_plus_0 = self.growth_factor (1.0)
        self.D_plus_i = self.D_plus (self.a_ini)
        self.mass_ini = 1.0
        self.mass_ini = self.mass (self.a_ini)
        self.expf_ini = 1.0
        self.expf_ini = self.expansion_function (self.a_ini)
    
    ## Calculates the exponent of an effective power law reproducing the
    #  initial slope of the growth factor
    def growth_index (self):
        e_equal = self.expansion_function (self.a_equal)
        v = [ \
             3.0+2.0*self.a_equal* \
             self.d_expansion_function_d_a (self.a_equal)/e_equal, \
             1.0-self.omega_m/(self.a_equal*(e_equal*self.a_equal)**2) \
            ]
        return 0.25*(-1.0-v[0]+np.sqrt ((1.0+v[0])**2+24.0*(1.0-v[1])))
    
    ## Returns the expansion function, that is the Hubble function divided by
    #  the Hubble constant, as a function of the scale factor \f$a\f$.
    #  For the cosmological standard model, the expansion function is
    #  \f[
    #    E(a) = \left(
    #      \Omega_ra^{-4}+\Omega_ma^{-3}+\Omega_ca^{-2}+\Omega_\Lambda
    #    \right)^{1/2}
    #  \f]
    def expansion_function (self, a):
        return np.sqrt (self.omega_r/a**4+
                        self.omega_m/a**3+
                        self.omega_l+
                        self.omega_c/a**2)
    
    ## Returns the derivative of the expansion function with respect to the
    #  scale factor \f$a\f$.
    #  For the cosmological standard model, the derivative of the expansion
    #  function is
    #  \f[
    #    \frac{dE}{da} = \frac{1}{2a^3E(a)}\left(
    #      -4\Omega_ra^{-2}-3\Omega_ma^{-1}-2\Omega_c
    #    \right)
    #  \f]
    def d_expansion_function_d_a (self, a):
        return 0.5/(a**3*self.expansion_function (a))* \
            (-4.0*self.omega_r/a**2-3.0*self.omega_m/a-2.0*self.omega_c)
    
    ## Returns the radiation-density parameter \f$\Omega_{r}\f$ as a function
    #  of the scale factor \f$a\f$
    #  \f[
    #    \Omega_r(a) = \frac{\Omega_{r0}}{a^4E^2(a)}
    #  \f]
    def omega_radiation (self, a):
        return self.omega_r/(a**4*self.expansion_function (a)**2)
    
    ## Returns the matter-density parameter \f$\Omega_{m}\f$ as a function
    #  of the scale factor \f$a\f$
    #  \f[
    #    \Omega_m(a) = \frac{\Omega_{m0}}{a^3E^2(a)}
    #  \f]
    def omega_matter (self, a):
        return self.omega_m/(a**3*self.expansion_function (a)**2)
    
    ## Returns the density parameter \f$\Omega_\Lambda\f$ of the cosmological
    #  constant as a function of the scale factor \f$a\f$
    #  \f[
    #    \Omega_\Lambda(a) = \frac{\Omega_{\Lambda0}}{E^2(a)}
    #  \f]
    def omega_lambda (self, a):
        return self.omega_l/(self.expansion_function (a)**2)
    
    ## Returns the cosmic time \f$t\f$ as a function of the scale factor \f$a\f$,
    #  in units of the Hubble time \f$H_0^{-1}\f$
    #  \f[
    #    t(a) = H_0^{-1}\int_0^a\frac{dx}{xE(x)}
    #  \f]
    def cosmic_time (self, a):
        return quad (lambda x: 1.0/(x*self.expansion_function (x)), 0, a)[0]
    
    ## Returns the radial metric distance between objects at redshifts \f$z_1\f$
    #  and \f$z_2\f$, in units of the Hubble radius \f$cH_0^{-1}\f$
    #  \f[
    #    r(z_1,z_2) = \frac{c}{H_0}\int_{a_1}^{a_2}\frac{dx}{x^2E(x)}\;,\quad
    #    a_i = (1+z_i)^{-1}
    #  \f]
    def radial_distance (self, z1, z2):
        a1 = 1.0/(1.0+z1)
        a2 = 1.0/(1.0+z2)
        return quad \
            (lambda x: 1.0/(x**2*self.expansion_function (x)), a1, a2)[0]
    
    ## Returns the angular metric distance \f$f_K(w)\f$ in units of the Hubble
    #  radius \f$cH_0^{-1}\f$
    #  \f[
    #    f_K(w) = \left\{
    #      \begin{array}{ll}
    #        \sinh w & \Omega_c < 0 \\ w & \Omega_c = 0 \\ \sin w & \Omega_c > 0
    #      \end{array}
    #    \right.
    #  \f]
    def f_K (self, w):
        if (abs (self.omega_c) <= 1.e-2):
            return w
        elif (self.omega_c < 0.0):
            return np.sinh (w)
        else:
            return np.sin (w)
        
    ## Returns the angular-diameter distance \f$D_A\f$ in units of the Hubble
    #  radius \f$cH_0^{-1}\f$
    #  \f[
    #    D_A(z_1,z_2) = \frac{1+z_2}{1+z_1}f_K(w)
    #  \f]
    def angular_distance (self, z1, z2):
        w = self.radial_distance (z1, z2)
        return (1.0+z2)/(1.0+z1)*self.f_K (w)
    
    ## Returns the luminosity distance \f$D_L\f$ in units of the Hubble
    #  radius \f$cH_0^{-1}\f$
    #  \f[
    #    D_L(z_1,z_2) = \frac{1+z_1}{1+z_2}f_K(w)
    #  \f]
    def luminosity_distance (self, z1, z2):
        w = self.radial_distance (z1, z2)
        return (1.0+z1)/(1.0+z2)*self.f_K (w)
    
    ## Implements the differential equation to be solved for the growth factor
    #  \f[
    #    \delta_k''+\left(\frac{3}{a}+\frac{E'(a)}{E(a)}\right)\delta_k'-
    #    \frac{3\Omega_m}{2a^2}\delta_k = 0
    #  \f]
    def growth_equation (self, delta, a):
        c = [ \
             1.5*self.omega_matter (a)/a**2, \
             3.0/a+self.d_expansion_function_d_a (a)/ \
             self.expansion_function (a) \
            ]
        return [delta[1], c[0]*delta[0]-c[1]*delta[1]]
    
    ## Returns the growth factor, which is the solution of the growth equation
    #  taken between scale factors \f$a_{eq}\f$ and \f$a\f$
    def growth_factor (self, a, n = 0):
        t = [self.a_equal, a]
        delta_initial = \
            [self.a_equal, self.n_growth*self.a_equal**(self.n_growth-1.0)]
        delta = odeint (self.growth_equation, delta_initial, t)
        return delta[1,n]
    
    ## Returns the growth factor, that is the solution of the growth equation
    #  taken between scale factors \f$a_{eq}\f$ and \f$a\f$, normalized to unity
    #  today at scale factor \f$a = 1\f$
    def D_plus (self, a):
        return self.growth_factor (a)/self.D_plus_0
    
    ## Returns the derivative of \f$D_+(a)\f$ with respect to the scale factor
    #  \f$a\f$
    def D_plus_prime (self, a):
        return self.growth_factor (a, 1)/self.D_plus_0
    
    ## Returns the effective KFT particle mass, normalized to unity at the initial
    #  scale factor
    #  \f[
    #    m(a) = a^2D_+\frac{d\ln D_+}{d\ln a}E = a^3D_+'E
    #  \f]
    def mass (self, a):
        return a**3*self.D_plus_prime (a)* \
            self.expansion_function (a)/self.mass_ini
    
    ## Returns the KFT time coordinate, that is \f$D_+/D_+^\mathrm{(i)}-1\f$
    def time (self, a):
        return self.D_plus (a)/self.D_plus_i-1
    
    ## Returns the Hamiltonian propagator
    #  \f[
    #    g_H(t',t) = \int_{t'}^t\frac{d\bar t}{m(\bar t)} =
    #      \int_{a'}^a\frac{d\bar t}{d\bar a}\frac{d\bar a}{m(\bar a)} =
    #      \int_{a'}^a\frac{dx}{x^3E(x)}
    #  \f]
    def g_H (self, a, a_prime = 1.0e-3):
        return quad \
          (lambda x: \
              (self.a_ini**3*self.expansion_function (self.a_ini))/ \
              (x**3*self.expansion_function (x)),
           a_prime, a)[0]/self.a_ini
