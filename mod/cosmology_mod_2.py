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
    def __init__ (self, c, a, \
        omega_radiation = 8.5e-5, \
        omega_matter = 0.3111, \
        omega_lambda = 0.6889, \
        hubble = 0.6766):
        # ,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15
        self.c        = c
        self.n_cheb   = len(self.c)
        # self.c1       = c1
        # self.c2       = c2
        # self.c3       = c3
        # self.c4       = c4
        # self.c5       = c5
        # self.c6       = c6
        # self.c7       = c7
        # self.c8       = c8
        # self.c9       = c9
        # self.c10      = c10
        # self.c11      = c11
        # self.c12      = c12
        # self.c13      = c13
        # self.c14      = c14
        # self.c15      = c15
        self.a        = a
        self.amin     = a[0]
        self.amax     = a[-1]
        self.a0       = a[0] # Too lazy to edit a0 to amin below

        self.omega_r  = omega_radiation
        self.omega_l  = omega_lambda
        self.omega_m  = omega_matter
        self.omega_b  = 0.04
        self.omega_c  = 1.0-self.omega_r-self.omega_m-self.omega_l
        self.hubble   = hubble

        # self.cheb_poly= np.polynomial.chebyshev.Chebyshev([self.c1,self.c2,self.c3,self.c4,self.c5], domain=[0,1])
        self.da       = (1-self.a[0])/self.a[0]
        self.x        = (self.a-self.a[0])/(self.a[-1]-self.a[0])
        self.dx_da    = 1.0 / (self.a[-1] - self.a[0])
        self.expf     = self.da*self.hubble/(1+self.da*self.x)**2/(self.cheb(self.c,2*self.x-1))

        P = self.cheb(self.c, 2*self.x-1)
        dP = self.cheb_derivative(self.c, 2*self.x-1)

        expf_raw = self.da * self.hubble / (1 + self.da * self.x)**2 / P
        self.d_expf_raw = (
            -2 * self.da**2 * self.hubble * self.dx_da / (1 + self.da * self.x)**3 / P
            - self.da * self.hubble * dP * self.dx_da / (1 + self.da * self.x)**2 / P**2
        )
        self.norm = expf_raw[-1]
        self.d_expf = (self.d_expf_raw * self.norm - expf_raw * self.d_expf_raw[-1]) / self.norm**2

        # self.rho_mean = 2.7594e11*self.hubble**2*self.omega_m
        self.a_equal  = self.omega_r/self.omega_m
        self.a_ini    = 1.0e-3
        self.n_growth = self.growth_index ()
        self.normalize ()
        # self.plot()
    
    def plot(self):
        np.savetxt("a_expf.txt",self.a)
        np.savetxt("expf_ori.txt", self.expansion_function(self.a))
    
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

    ## Using a_equal
    # def growth_index (self):
    #     e_equal = self.expansion_function (self.a_equal)
    #     # e_equal = 152234.18472252088 # Hard-coded
    #     v = [ \
    #          3.0+2.0*self.a_equal* \
    #          self.d_expansion_function_d_a (self.a_equal)/e_equal, \
    #          1.0-self.omega_m/(self.a_equal*(e_equal*self.a_equal)**2) \
    #         ]
    #     # return 0.7499922269714218 # Hard-coded
    #     return 0.25*(-1.0-v[0]+np.sqrt ((1.0+v[0])**2+24.0*(1.0-v[1])))

    # def growth_factor (self, a, n = 0):
    #     t = [self.a_equal, a]
    #     delta_initial = \
    #         [self.a_equal, self.n_growth*self.a_equal**(self.n_growth-1.0)]
    #     delta = odeint (self.growth_equation, delta_initial, t)
    #     return delta[1,n]
    
    ## Using a(z=6)
    def growth_index (self):
        e0 = self.expansion_function (self.a0)
        v = [ \
             3.0+2.0*self.a0* \
             self.d_expansion_function_d_a (self.a0)/e0, \
             1.0-self.omega_m/(self.a0*(e0*self.a0)**2) \
            ]
        return 0.25*(-1.0-v[0]+np.sqrt ((1.0+v[0])**2+24.0*(1.0-v[1])))
    
    def growth_factor (self, a, n = 0):
        t = [self.a0, a]
        delta_initial = \
            [self.a0, self.n_growth*self.a0**(self.n_growth-1.0)]
        delta = odeint (self.growth_equation, delta_initial, t)
        return delta[1,n]
    

    ## LCDM expansion function
    # def expansion_function (self, a):
    #     return np.sqrt (self.omega_r/a**4+
    #                     self.omega_m/a**3+
    #                     self.omega_l+
    #                     self.omega_c/a**2)
    
    # def d_expansion_function_d_a (self, a):
    #     return 0.5/(a**3*self.expansion_function (a))* \
    #         (-4.0*self.omega_r/a**2-3.0*self.omega_m/a-2.0*self.omega_c)
    
    ## Model-independent expansion function
    # def expansion_function (self, a):
    #     # cheb_poly = np.polynomial.chebyshev.Chebyshev([self.c1,self.c2,self.c3,self.c4,self.c5,self.c6,self.c7,self.c8,self.c9,self.c10,self.c11,self.c12,self.c13,self.c14,self.c15], domain=[0,1])
    #     cheb_poly = np.polynomial.chebyshev.Chebyshev([self.c1,self.c2,self.c3,self.c4,self.c5], domain=[0,1])
    #     da = (1-self.amin)/self.amin
    #     x = (a-self.amin)/(self.amax-self.amin)
    #     expf = da*self.hubble/(1+da*x)**2/(cheb_poly(x))
    #     # expf = expf/expf[-1]
    #     expf = expf/self.expf[-1]
    #     # expf = self.hubble*cheb_poly(a)
    #     return expf

    def cheb(self,c,a):
            print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ a = ",a)
            x = np.arccos(a)
            result = 0
            for i in range(len(c)) :    
                if i ==0 :
                    result += (1/np.pi)**(1/2)*c[i]
                else :
                    result += (2/np.pi)**(1/2)*c[i]*np.cos(i*x)
            return result

    def cheb_derivative(self,c,x) :
        n=len(c)
        x=np.atleast_1d(x)
        q=np.zeros((n,x.shape[0]))
        t=np.zeros((n,x.shape[0]))
        t[0]=1
        t[1]=x
        for i in range(1,n-1) :
            t[i+1]=2*x*t[i]-t[i-1]
        for i in range(n):
            if i ==0 :
                t[i] *= (1/np.pi)**(1/2)
            else :
                t[i] *= (2/np.pi)**(1/2)
        q[0]=0
        q[1]=(2/np.pi)**(1/2)
        for i in range(1,n-1) :
            q[i+1]=2*t[i] + 2*x*q[i] - q[i-1]
        result = 0
        for i in range(n) :
            result += c[i]*q[i]
        if result.shape[0] == 1 :
            result = result[0]
        return result
    
    def expansion_function (self, a):
        x = (a-self.amin)/(self.amax-self.amin)         
        expf = self.da*self.hubble/(1+self.da*x)**2/(self.cheb(self.c,2*x-1))
        # expf = expf/expf[-1]
        expf = expf/self.expf[-1]
        # expf = self.hubble*cheb_poly(a)
        return expf

    # Numerical derivative
    # def d_expansion_function_d_a (self, a):
    #     # cheb_poly = np.polynomial.chebyshev.Chebyshev([self.c1,self.c2,self.c3,self.c4,self.c5,self.c6,self.c7,self.c8,self.c9,self.c10,self.c11,self.c12,self.c13,self.c14,self.c15], domain=[0,1])
    #     cheb_poly = np.polynomial.chebyshev.Chebyshev([self.c1,self.c2,self.c3], domain=[0,1])
    #     da = (1-a[0])/a[0]
    #     x = (a-a[0])/(a[-1]-a[0])
    #     expf = da*self.hubble/(1+da*x)**2/(cheb_poly(x))
    #     expf = expf/expf[-1]
    #     d_expf = np.gradient(expf,a)
    #     # cheb_deriv = cheb_poly.deriv()
    #     # d_expf = self.hubble*cheb_deriv(a)
    #     return d_expf
    
    # Analytical derivative
    # def d_expansion_function_d_a (self, a):
    #     cheb_poly = np.polynomial.chebyshev.Chebyshev([self.c1,self.c2,self.c3,self.c4,self.c5], domain=[0,1])
    #     da = (1-self.amin)/self.amin
    #     x = (a-self.amin)/(self.amax-self.amin)
    #     dx_da = 1.0 / (self.amax - self.amin)

    #     P = cheb_poly(x)
    #     dP = cheb_poly.deriv()(x)

    #     expf_raw = da * self.hubble / (1 + da * x)**2 / P
    #     d_expf_raw = (
    #         -2 * da**2 * self.hubble * dx_da / (1 + da * x)**3 / P
    #         - da * self.hubble * dP * dx_da / (1 + da * x)**2 / P**2
    #     )
    #     norm = self.expf[-1]
    #     d_expf = (d_expf_raw * norm - expf_raw * self.d_expf_raw[-1]) / norm**2
    #     return d_expf
    def d_expansion_function_d_a (self, a):
        # da = (1-self.amin)/self.amin
        x = (a-self.amin)/(self.amax-self.amin)
        # dx_da = 1.0 / (self.amax - self.amin)

        P = self.cheb(self.c, 2*x-1)
        dP = self.cheb_derivative(self.c, 2*x-1)

        expf_raw = self.da * self.hubble / (1 + self.da * x)**2 / P
        d_expf_raw = (
            -2 * self.da**2 * self.hubble * self.dx_da / (1 + self.da * x)**3 / P
            - self.da * self.hubble * dP * self.dx_da / (1 + self.da * x)**2 / P**2
        )
        norm = self.expf[-1]
        d_expf = (d_expf_raw * norm - expf_raw * self.d_expf_raw[-1]) / norm**2
        return d_expf
    
    ## Implements the differential equation to be solved for the growth factor
    def growth_equation (self, delta, a):
        c = [ \
             1.5*self.omega_matter (a)/a**2, \
             3.0/a+self.d_expansion_function_d_a (a)/ \
             self.expansion_function (a) \
            ]
        return [delta[1], c[0]*delta[0]-c[1]*delta[1]]


    ## Returns the growth factor, that is the solution of the growth equation
    #  taken between scale factors \f$a_{eq}\f$ and \f$a\f$, normalized to unity
    #  today at scale factor \f$a = 1\f$
    def D_plus (self, a):
        return self.growth_factor (a)/self.D_plus_0
    
    def omega_matter (self, a):
        return self.omega_m/(a**3*self.expansion_function (a)**2)

##########################################################################################################################3

    # ## Returns the radiation-density parameter \f$\Omega_{r}\f$ as a function
    # #  of the scale factor \f$a\f$
    # #  \f[
    # #    \Omega_r(a) = \frac{\Omega_{r0}}{a^4E^2(a)}
    # #  \f]
    # def omega_radiation (self, a):
    #     return self.omega_r/(a**4*self.expansion_function (a)**2)
    

    
    # ## Returns the density parameter \f$\Omega_\Lambda\f$ of the cosmological
    # #  constant as a function of the scale factor \f$a\f$
    # #  \f[
    # #    \Omega_\Lambda(a) = \frac{\Omega_{\Lambda0}}{E^2(a)}
    # #  \f]
    # def omega_lambda (self, a):
    #     return self.omega_l/(self.expansion_function (a)**2)
    
    # ## Returns the cosmic time \f$t\f$ as a function of the scale factor \f$a\f$,
    # #  in units of the Hubble time \f$H_0^{-1}\f$
    # #  \f[
    # #    t(a) = H_0^{-1}\int_0^a\frac{dx}{xE(x)}
    # #  \f]
    # def cosmic_time (self, a):
    #     return quad (lambda x: 1.0/(x*self.expansion_function (x)), 0, a)[0]
    
    # ## Returns the radial metric distance between objects at redshifts \f$z_1\f$
    # #  and \f$z_2\f$, in units of the Hubble radius \f$cH_0^{-1}\f$
    # #  \f[
    # #    r(z_1,z_2) = \frac{c}{H_0}\int_{a_1}^{a_2}\frac{dx}{x^2E(x)}\;,\quad
    # #    a_i = (1+z_i)^{-1}
    # #  \f]
    # def radial_distance (self, z1, z2):
    #     a1 = 1.0/(1.0+z1)
    #     a2 = 1.0/(1.0+z2)
    #     return quad \
    #         (lambda x: 1.0/(x**2*self.expansion_function (x)), a1, a2)[0]
    
    # ## Returns the angular metric distance \f$f_K(w)\f$ in units of the Hubble
    # #  radius \f$cH_0^{-1}\f$
    # #  \f[
    # #    f_K(w) = \left\{
    # #      \begin{array}{ll}
    # #        \sinh w & \Omega_c < 0 \\ w & \Omega_c = 0 \\ \sin w & \Omega_c > 0
    # #      \end{array}
    # #    \right.
    # #  \f]
    # def f_K (self, w):
    #     if (abs (self.omega_c) <= 1.e-2):
    #         return w
    #     elif (self.omega_c < 0.0):
    #         return np.sinh (w)
    #     else:
    #         return np.sin (w)
        
    # ## Returns the angular-diameter distance \f$D_A\f$ in units of the Hubble
    # #  radius \f$cH_0^{-1}\f$
    # #  \f[
    # #    D_A(z_1,z_2) = \frac{1+z_2}{1+z_1}f_K(w)
    # #  \f]
    # def angular_distance (self, z1, z2):
    #     w = self.radial_distance (z1, z2)
    #     return (1.0+z2)/(1.0+z1)*self.f_K (w)
    
    # ## Returns the luminosity distance \f$D_L\f$ in units of the Hubble
    # #  radius \f$cH_0^{-1}\f$
    # #  \f[
    # #    D_L(z_1,z_2) = \frac{1+z_1}{1+z_2}f_K(w)
    # #  \f]
    # def luminosity_distance (self, z1, z2):
    #     w = self.radial_distance (z1, z2)
    #     return (1.0+z1)/(1.0+z2)*self.f_K (w)
    
    
    # ## Returns the derivative of \f$D_+(a)\f$ with respect to the scale factor
    # #  \f$a\f$
    def D_plus_prime (self, a):
        return self.growth_factor (a, 1)/self.D_plus_0
    
    # ## Returns the effective KFT particle mass, normalized to unity at the initial
    # #  scale factor
    # #  \f[
    # #    m(a) = a^2D_+\frac{d\ln D_+}{d\ln a}E = a^3D_+'E
    # #  \f]
    def mass (self, a):
        return a**3*self.D_plus_prime (a)* \
            self.expansion_function (a)/self.mass_ini
    
    # ## Returns the KFT time coordinate, that is \f$D_+/D_+^\mathrm{(i)}-1\f$
    # def time (self, a):
    #     return self.D_plus (a)/self.D_plus_i-1
    
    # ## Returns the Hamiltonian propagator
    # #  \f[
    # #    g_H(t',t) = \int_{t'}^t\frac{d\bar t}{m(\bar t)} =
    # #      \int_{a'}^a\frac{d\bar t}{d\bar a}\frac{d\bar a}{m(\bar a)} =
    # #      \int_{a'}^a\frac{dx}{x^3E(x)}
    # #  \f]
    # def g_H (self, a, a_prime = 1.0e-3):
    #     return quad \
    #       (lambda x: \
    #           (self.a_ini**3*self.expansion_function (self.a_ini))/ \
    #           (x**3*self.expansion_function (x)),
    #        a_prime, a)[0]/self.a_ini

    # def cheb_derivative(self, c, a):
    #     x = np.arccos(a)
    #     dx_da = -1 / np.sqrt(1 - a**2)  # derivative of arccos(a)

    #     result = 0
    #     for i in range(1, len(c)):  # derivative of constant term (i=0) is zero
    #         coeff = (2 / np.pi)**0.5 * c[i]
    #         d_cosix_dx = -i * np.sin(i * x)  # derivative of cos(i * x)
    #         result += coeff * d_cosix_dx
    #     return result * dx_da
    
