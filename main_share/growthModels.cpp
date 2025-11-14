#include <iostream>
#include <astro/utilities/utilities.h>
#include <astro/utilities/deSolver.h>

#include "growthModels.h"
#include <astro/cosmology/standardGrowthFactor.h>
#include <astro/cosmology/cosmologyBase.h>

growthModels::growthModels
  (expansionModel * exp_model_in, double a_min_in,
   double omega_m0_in):
exp_model (exp_model_in), a_min (a_min_in), omega_m0 (omega_m0_in)
{

  growth_tab.fill
    ([&] (double a, double omega_m0)
    {
      return growth(a, omega_m0);

    }, astro::LOG_SPACING, astro::LINEAR_SPACING, 1028, 1028,
    a_min, 1.0, 0.0, 1.0);

    std::cout<<"finished table"<<std::endl;

  // deltaDplus_tab.fill
  //   ([&] (double a, double omega_m0)
  //   {
  //     return deltaDplus(a, omega_m0);
  //
  //   }, astro::LOG_SPACING, astro::LINEAR_SPACING, 64, 1,
  //   a_min, 1.0, 0.1, 1.0);

    std::cout<<"finished table 2"<<std::endl;

}

double growthModels::operator () (double a, double k)
{ return a*reducedGrowth (a,omega_m0); }

double growthModels::Dplus (double a)
{ return growth_tab(a,omega_m0)/growth_tab(1.0,omega_m0); }

double growthModels::Dplus_marg (double a, double mean, double sigma)
{
  astro::integrator integral ([&] (double y)
  { return
      growth_tab(a,y)/growth_tab(1.0,y)*gaussFilter(y,mean,sigma); });

  astro::integrator norm ([&] (double y)
  { return
      gaussFilter(y,mean,sigma); });

  return integral(0.0,0.7)/norm(0.0,0.7); }

double growthModels::Dplus_marg_std (double a, double mean, double sigma)
{
  astro::integrator integral ([&] (double y)
  { return
      pow(growth_tab(a,y)/growth_tab(1.0,y)-Dplus_marg(a,mean,sigma),2.0)*gaussFilter(y,mean,sigma); },1.0e-3);

  astro::integrator norm ([&] (double y)
  { return
      gaussFilter(y,mean,sigma); });
  return integral(0.0,0.7)/norm(0.0,0.7); }

double growthModels::deltaDplus_marg (double a, double mean, double sigma)
{
  astro::integrator integral ([&] (double y)
  { return
      deltaDplus(a,y)*gaussFilter(y,mean,sigma); });

  astro::integrator norm ([&] (double y)
  { return
      gaussFilter(y,mean,sigma); });

  return integral(0.0,0.7)/norm(0.0,0.7); }

double growthModels::gaussFilter (double x, double mu, double sigma){

  return 1./sigma/sqrt(2.*M_PI)*exp(-1./2.*gsl_pow_2((x-mu)/sigma));
}

double growthModels::dlnDplus_dlna (double a)
{
  const double epsilon = 1.0e-3;
  double a_plus = std::min(1.0,(1.0+epsilon)*a);;
  double a_minus = (1.0-epsilon)*a;
  double Dplus_plus = growth_tab(a_plus, omega_m0);
  double Dplus_minus = growth_tab(a_minus, omega_m0);
  double delta_Dplus = (Dplus_plus-Dplus_minus);
  return a/growth_tab(a, omega_m0)*delta_Dplus/(a_plus-a_minus);
}

double growthModels::dDplus_dOmega_m02 (double a)
{
  const double epsilon = 1.0e-3;
  double Om_plus = (1.0+epsilon)*omega_m0;
  double Om_minus = (1.0-epsilon)*omega_m0;
  double Dplus_plus = growth_tab(a, Om_plus)/growth_tab(1.0,Om_plus);
  double Dplus_minus = growth_tab(a, Om_minus)/growth_tab(1.0,Om_minus);
  return (Dplus_plus-Dplus_minus)/(Om_plus-Om_minus);
}

double growthModels::dDplus_dOmega_m0 (double a)
{
  astro::integrator int1 ([&] (double y)
    { return
        1.0/exp_model->expansionFunction (y)/pow(y, 2.0)
        *gsl_pow_2 (growth_tab(y,omega_m0)); });

  astro::integrator int2 ([&] (double y)
    { double e  = exp_model->expansionFunction (a_min);
      double de = exp_model->dexpansionFunction_da (a_min);
      double dn_dOm = 3./4./(index(a_min,omega_m0)+1.+de/e*a_min/2.)
                      /pow(a_min,3.)/pow(e,2.);
      double k0 = exp_model->expansionFunction (a_min)*pow(a_min, 2.0)*dn_dOm;
      return
        1.0/exp_model->expansionFunction (y)/pow(y, 3.0)/
        gsl_pow_2 (growth_tab(y,omega_m0))*(3./2.*int1(a_min,y)+k0); });

  return Dplus(a)*int2(1.0, a);

  // astro::deSolver solve ([&] (double a, const double y[], double dy[])
  // {
  //   double e  = exp_model->expansionFunction (a);
  //   double de = exp_model->dexpansionFunction_da (a);
  //   double Omega = omega_m0/pow (a*e, 2)/a;
  //   dy[0] = y[1];
  //   dy[1] = -(3.0/a+de/e)*y[1]+3.0*Omega/2.0/pow (a, 2)*y[0]+3./2./pow(a,5.)/pow(e,2.)*growth(a,omega_m0);
  //   return 0;
  // }, 2, 1.0e-4);
  //
  // double e  = exp_model->expansionFunction (a_min);
  // double de = exp_model->dexpansionFunction_da (a_min);
  // double dn_dOm = 3./4./(index(a_min,omega_m0)+1.+de/e*a_min/2.)/pow(a_min,3.)/pow(e,2.);
  // double g_init[2] = {0.0,dn_dOm/a_min};
  // std::vector<double> g_max = solve (g_init, a_min, 1.0);
  // std::vector<double> g = solve (g_init, a_min, a);
  // return g[0]/growth(1.0,omega_m0)-Dplus(a)/growth(1.0,omega_m0)*g_max[0];
}

double growthModels::index (double a, double Omega_m0){

  double e = 3.+2.*exp_model->dexpansionFunction_da (a)
             /exp_model->expansionFunction(a)*a;
  double w = 1.-Omega_m0/pow(exp_model->expansionFunction(a),2.0)/pow(a,3.0);

  return 1./4.*(-1.-e+sqrt(pow(1.+e,2.)+24.*(1.-w)));
}

double growthModels::growth (double a,double Omega_m0)
{
  astro::deSolver solve ([&] (double a, const double y[], double dy[])
  {
    double e  = exp_model->expansionFunction (a);
    double de = exp_model->dexpansionFunction_da (a);
    double Omega = Omega_m0/pow (a*e, 2)/a;
    dy[0] = y[1];
    dy[1] = -(3.0/a+de/e)*y[1]+3.0*Omega/2.0/pow (a, 2)*y[0];
    return 0;
  }, 2, 1.0e-4);

  double g_init[2] = {1.0, index(a_min,Omega_m0)/a_min};
  std::vector<double> g_max = solve (g_init, a_min, 1.0);
  std::vector<double> g = solve (g_init, a_min, a);
  return g[0];

}

double growthModels::reducedGrowth (double a, double Omega_m0)
{
  astro::deSolver solve ([&] (double a, const double y[], double dy[])
  {
    double e  = exp_model->expansionFunction (a);
    double de = exp_model->dexpansionFunction_da (a);
    double Omega = Omega_m0/pow (a*e, 2)/a;
    dy[0] = y[1];
    dy[1] = -(5.0/a+de/e)*y[1]-(3.0*(1.0-Omega/2.0)/a+de/e)*y[0]/a;
    return 0;
  }, 2, 1.0e-4);

  double g_init[2] = {1.0, (index(a_min,Omega_m0)-1.0)/a_min};
  std::vector<double> g_max = solve (g_init, a_min, 1.0);
  std::vector<double> g = solve (g_init, a_min, a);
  return g[0]/g_max[0];
}

double growthModels::gravitationalConstant (double a, double k)
{ return 1.0; }

double growthModels::dDplus_dE (double x, double a, double Omega_m0)
{
  astro::integrator Gamma ([&] (double y)
    { return
        1.0/exp_model->expansionFunction (y)/y/
        gsl_pow_2 (y*growth_tab(y, Omega_m0)); });

  double E       = exp_model->expansionFunction (x);
  double Omega_m = Omega_m0/x/gsl_pow_2 (x*E);
  // double gamma   = exp_model->gamma (x, omega_m0, NULL);
  // std::cout<<gamma<<std::endl;
  // gamma = 6.0/11.0;
  // double f_of_x  = x*gsl_pow_2 (Dplus (x))*
  //   Omega_m*(pow (Omega_m, 2.0*gamma-1.0)-1.5);
  //double g_of_x  = pow (Omega_m, gamma)/x/E;
  double f_of_x  = x*gsl_pow_2 (growth_tab(x,Omega_m0))*
    (pow (dlnDplus_dlna(x), 2.0)-Omega_m*1.5);
  double g_of_x  = dlnDplus_dlna(x)/x/E;

  return
    growth_tab(a, Omega_m0)*(
      ((x <= a) ? (f_of_x*Gamma (x, a)
      -g_of_x) : 0.0)
    );
}

double growthModels::deltaE (double a)
{ return exp_model->deltaE (a); }

double growthModels::deltaDplus (double a, double Omega_m0)
{
  astro::integrator error ([&] (double x)
    { return 1.0/growth_tab(1.0, Omega_m0)*(dDplus_dE (x, a, Omega_m0))
      *deltaE (x); },1.0e-3);
  astro::integrator error2 ([&] (double x)
    { return 1.0/growth_tab(1.0, Omega_m0)*(-growth_tab(a,Omega_m0)/growth_tab(1.0,Omega_m0)
      *dDplus_dE(x,1.0, Omega_m0))*deltaE (x); },1.0e-3);
  return error (a_min, a)+error2(a_min,1.0);
}

void growthModels::print (std::string filename, int n)
{
  astro::cosmologyBase cos (omega_m0,1.0-omega_m0,0.7,0.0);
  cos.setDarkUniverse ();
  astro::standardGrowthFactor growth_lcdm (&cos);

  std::ofstream output (filename);
  for (int i = 0; i < n; i++)
  {
    double a = astro::x_logarithmic (i, n, a_min, 1.0);
    output
      << a
      << "  " << Dplus (a)
      << "  " << deltaDplus (a,omega_m0)
      << "  " << growth_lcdm (a)/growth_lcdm(1.0)*a
      << "  " << dDplus_dOmega_m0(a)
      << "  " << Dplus_marg (a,omega_m0,0.05)
      << "  " << deltaDplus_marg (a,omega_m0,0.05)
      << "  " << Dplus_marg_std(a,omega_m0,0.05)
      << std::endl;
  }
  output.close ();
}
