import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

class model_independent:
    def __init__(self,c,a,h0,om,n) :
        self.a = a
        self.amin = self.a[0]
        self.amax = self.a[-1]
        self.da = (1-a[0])/a[0]
        self.x = (a-a[0])/(a[-1]-a[0])
        self.h0=h0
        self.om = om
        self.c = c
        self.n = n
        self.norm = self._compute_norm()
    
    def _compute_norm(self):
        expf = self.da/self.amin**(self.n/2)/(1+self.x*self.da)**(2+self.n/2)/self.chebs(self.c,2*self.x-1)/self.h0
        return expf[-1]

    # This is a method to call Chebyshev polynomials of order given by len(c)
    def chebs (self,c,a) :
        if np.ndim(a) == 0 : # Needed for odeint, for some reason it asigns values outside boundaries
            if a>1 :
                a=1.0
            elif  a<-1 :
                a=-1.0
        x = np.arccos(a)
        result = 0
        for i in range(len(c)) :
            if i ==0 :
                result += (1/np.pi)**(1/2)*c[i]
            else :
                result += (2/np.pi)**(1/2)*c[i]*np.cos(i*x)
        return result
    
    # This is also a method to call Chebyshev polynomials of order given by len(c), similar to above
    # I am simply testing whether the recursive method works, it does.
    def cheb_recursive(c,x) :
        n = len(c)
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
        result = 0
        for i in range(len(c)) :
            result += c[i]*t[i]
        return result
    
    # This is a method to call just the n-th order of the Chebyshev polynomials,
    # not the be confused with the previous methods that call the entire polynomials
    def cheb(self,n,x) :
        n = n+1 # For n-th polynomial, we need n+1 eements because the polynomial starts from 0
        t=np.zeros((n,x.shape[0]))
        # t[0]=(1/np.pi)**(1/2)
        # t[1]=(2/np.pi)**(1/2)*x
        if n == 1 :
            t[0] = (1/np.pi)**(1/2)
            return t[0]
        t[0]=1
        t[1]=x
        for i in range(1,n-1) :
            t[i+1]=2*x*t[i]-t[i-1]
        for i in range(n):
            if i ==0 :
                t[i] *= (1/np.pi)**(1/2)
            else :
                t[i] *= (2/np.pi)**(1/2)
        return t[n-1]
    
    # Method to call the derivative of the Chebyshev polynomials of order given by len(c)
    def cheb_derivative(self,c,x) :
        n=len(c)
        x=np.atleast_1d(x) # Needed for odeint
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
        if result.shape[0] == 1 : # Needed for odeint
            result = result[0]
        return result
    
    # This gives the LCDM expansion function with cosmology used to define the class
    def lcdm_expf (self) :
        omega_m = self.om
        omega_l = 1-omega_m
        expf = np.sqrt (omega_m/self.a**3 + omega_l)
        return expf

    # This gives the derivative of the LCDM expansion function
    def d_expansion_function_d_a_lcdm (self, a):
        return 0.5/(a**3*self.lcdm_expf ())* \
            (-3.0*self.om/a)
        # return 0.5/(a**3*self.lcdm_expf (a))* \
        #     (-4.0*self.omega_r/a**2-3.0*self.omega_m/a-2.0*self.omega_c)
    
    # This gives the model-independent expansion function (Eq. 4.3.1 in my thesis)
    def model_expf (self, a=None, c=None, n=None) :
        if n is None :
            n = self.n
        if c is None :
            c = self.c
        if a is not None : # Needed for odeint
            x = (a-self.amin)/(self.amax-self.amin)
        else :
            x = self.x
        expf = self.da/self.amin**(n/2)/(1+x*self.da)**(2+n/2)/self.chebs(c,2*x-1)/self.h0
        expf = expf/self.norm
        return expf
    
    # This gives the derivative of the model-independent expansion function (derivative of Eq. 4.3.1 in my thesis)
    def d_expansion_function_d_a (self,a=None,c=None,n=None):
        if n is None :
            n = self.n
        if c is None :
            c = self.c
        if a is not None : # Needed for odeint
            x = (a-self.amin)/(self.amax-self.amin)
        else :
            x = self.x
        A = (1+x*self.da)
        prefactor = -self.da/self.h0/self.amin**(n/2)/(self.a[-1]-self.a[0])/A**(2+n/2)/self.chebs(c,2*x-1)
        term1 = (2+n/2)*self.da/A
        term2 = 2*self.cheb_derivative(c,2*x-1)/self.chebs(c,2*x-1)
        d_expf = prefactor*(term1 + term2)/self.norm
        return d_expf
    
    # This calculates the comoving distance with the model-independent expansion function
    def fast_comoving_distance(self, a_vals=None, expf_vals=None, h0=None):
        if a_vals is None:
            a_vals = self.a
        if expf_vals is None:
            expf_vals = self.model_expf()
        if h0 is None:
            h0 = self.h0
        c_light = 299792.458  # km/s
        integrand = 1 / (a_vals**2 * expf_vals)

        # Flip arrays to integrate from low to high a
        a_flip = a_vals[::-1]
        integrand_flip = integrand[::-1]

        cum_int = cumulative_trapezoid(integrand_flip, a_flip, initial=0)

        # Reverse result to match original order of a_vals, multiply by -1 to get positive result
        d_c = -c_light / (h0 * 100) * cum_int[::-1]
        return d_c
    
    # This gives the expansion function using a PCA model, where the Chebyshev polynomials define the deviation from LCDM
    # (USE WITH CAUTION, NOT TESTED!!)
    def model_pca(self,c=None) :
        if c is None :
            c = self.c
        expf = self.lcdm_expf() + self.chebs(c,2*self.x-1)
        return expf
    
    # This fits the model-independent expansion function (Eq. 4.3.1) into the LCDM expansion function using gradient descent method
    # It then updates the coefficients 'c' in the class to the fitted values
    def fit_model_to_lcdm (self,model_init,learning_rate = 0.1,n=10) :
        y_hat = model_init
        y = self.lcdm_expf()
        loss = np.sum(1/2*(y-y_hat)**2)
        print(loss)
        grad = np.zeros(len(self.c))
        grad_y = np.zeros((len(self.c),self.x.shape[0]))
        
        for i in range(0,n+1) :
            for j in range(0,len(self.c)) :
                grad_y[j] = -self.da/self.h0 * self.cheb(j,2*self.x-1)/(1+self.da*self.x)**2/self.chebs(self.c,2*self.x-1)**2
                grad[j] = np.sum((y_hat-y)*grad_y[j])

            self.c = self.c - learning_rate*grad
            self.norm = self._compute_norm()
            y_hat = self.model_expf()
            loss = np.sum(1/2*(y-y_hat)**2)
            if i == 0 or i%1000==0:
                print(f"iteration {i} ",loss)

    # This fits the model-independent expansion function given by the 'PCA' method into the LCDM expansion function using gradient descent method
    # It then updates the coefficients 'c' in the class to the fitted values
    # (USE WITH CAUTION, NOT TESTED!!)
    def fit_model_to_lcdm_pca (self,model_init,learning_rate = 0.1,n=10) :
        y = self.lcdm_expf()
        y_hat = model_init
        loss = np.sum(1/2*(y-y_hat)**2)
        print(loss)
        grad = np.zeros(len(self.c))
        grad_y = np.zeros((len(self.c),self.x.shape[0]))
        
        for i in range(0,n+1) :
            for j in range(0,len(self.c)) :
                grad_y[j] = self.cheb(j,2*self.x-1)
                grad[j] = np.sum((y_hat-y)*grad_y[j])

            self.c = self.c - learning_rate*grad
            self.norm = self._compute_norm()
            y_hat = self.model_pca()
            loss = np.sum(1/2*(y-y_hat)**2)
            if i == 0 or i%1000==0:
                print(f"iteration {i} ",loss)
                print(f"c : {self.c} \n")

# Class to calculate the growth factor model independently
# Adapted from Matthias' Cosmology class
class Growth_Factor :
    def __init__(self,c,a,h0,om,n):
        self.c = c
        self.a = a
        self.amin = self.a[0]
        self.amax = self.a[-1]
        self.n = n
        self.model = model_independent(self.c,self.a,h0,om,self.n)
        self.n_growth = self.growth_index ()
        self.normalize()
    
    # This gives Eq. 6.5 in Alex' thesis
    def growth_index (self):
        e0 = self.model.model_expf (a=self.amin)
        v = [ \
            3.0+2.0*self.amin* \
            self.model.d_expansion_function_d_a (a=self.amin)/e0, \
            1.0-self.model.om/(self.amin*(e0*self.amin)**2) \
            ]
        return 0.25*(-1.0-v[0]+np.sqrt ((1.0+v[0])**2+24.0*(1.0-v[1])))
    
    def normalize (self):
        self.D_plus_0 = 1.0
        self.D_plus_0 = self.growth_factor (1.0)

    def omega_matter (self, a):
        return self.model.om/(a**3*self.model.model_expf (a=a)**2)
    
    def growth_factor (self, a, n = 0):
        t = [self.amin, a]
        delta_initial = \
            [self.amin, self.n_growth*self.amin**(self.n_growth-1.0)]
        delta = odeint (self.growth_equation, delta_initial, t)
        return delta[1,n]
    
    def growth_equation (self, delta, a):
        c = [ \
            1.5*self.omega_matter (a)/a**2, \
            3.0/a+self.model.d_expansion_function_d_a (a=a)/ \
            self.model.model_expf (a=a) \
        ]
        return [delta[1], c[0]*delta[0]-c[1]*delta[1]]
    
    def D_plus (self, a):
        return self.growth_factor (a)/self.D_plus_0

# Class to plot graphs
# (NOT WELL-TESTED!)
class Report :
    def __init__(self,om,h0,n):
        self.om = om
        self.h0 = h0
        self.n = n

    def combine(self,c,theory,pm,cell, comment="") :
        k = np.loadtxt("k_hme.txt")
        a = np.loadtxt("a_bg.txt")
        ells = np.loadtxt("ellsme.txt")
        pm_lcdm = np.loadtxt("pm_lin_60.txt")[0]
        cell_lcdm = np.loadtxt("c_ell_lin_60.txt")
        self.model = model_independent(c,a,self.h0,self.om,self.n)
        model = self.model.model_expf(c=c)
        lcdm = self.model.lcdm_expf()
        plt.figure(figsize=(15,5))
        plt.subplot(1,4,1)
        plt.plot(a,lcdm,label="lcdm")
        plt.plot(a,model,label="model")
        # plt.xlim(0.5,1)
        plt.ylabel("E(a)")
        plt.xlabel("a")
        plt.title(comment)
        plt.legend()

        dt60 = np.loadtxt("data-theory_data60_theory60.txt")
        diff_best_like = np.abs(theory-dt60[:,0])
        diff_lcdm_60 = np.abs(dt60[:,1]-dt60[:,0])
        n=40
        print("difference from chain", comment, ":", np.sum(diff_best_like[n:]))
        print("difference from lcdm 60", comment, ":" , np.sum(diff_lcdm_60[n:]))
        plt.subplot(1,4,2)
        plt.plot(dt60[:,1],label="theory_lcdm",marker = '.')
        plt.plot(theory, label = "theory_cheb",marker='.')
        plt.plot(dt60[:,0], label = "data",marker = '.')
        plt.title("COSEBIs (En)")
        plt.legend();

        plt.subplot(1,4,3)
        plt.loglog(k,pm_lcdm,label="lcdm")
        plt.loglog(k,pm,label="cheb")
        plt.legend();
        plt.xlabel('k')
        plt.title("matter power spectrum")

        plt.subplot(1,4,4)
        plt.loglog(ells,cell_lcdm,label="lcdm")
        plt.loglog(ells,cell,label="cheb")
        plt.legend();
        plt.xlabel('ell')
        plt.title("Shear power spectrum bin 5,5")