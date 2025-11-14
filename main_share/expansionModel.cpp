#include "expansionModel.h"

#include <iostream>
#include <astro/utilities/utilities.h>
#include <astro/utilities/nlModeller.h>
#include <astro/utilities/integrator.h>
#include <astro/cosmology/cosmologyBase.h>

expansionModel::expansionModel
  (astro::functionSet * basis_in, std::vector<double>& c_in,
   std::vector<double>& dc_in, double a_min_in, double exponent_in):
basis (basis_in), c (c_in), dc (dc_in), a_min (a_min_in), exponent (exponent_in)
{ init_model (); }

void expansionModel::init_model ()
{
  n_basis = c.size ();
  delta_a = 1.0/a_min-1.0;
  norm = 1.0;
  norm = 1.0/expansionFunction (1.0);

  astro::cosmologyBase cos_model;
  Omega_r0 = 0.0; // cos_model.omegaRadiation ();
  fitModel ();
}

double expansionModel::mf (double a, std::vector<double> p)
{ return p[0]+p[1]*gsl_pow_3 (a); }

std::vector<double> expansionModel::dm
  (double a, std::vector<double> p)
{
  std::vector<double> df;
  df.push_back (1.0);
  df.push_back (gsl_pow_3 (a));
  return df;
}

void expansionModel::fitModel ()
{
  astro::data_model data;
  data.n = 32;
  data.p = 2;
  data.x = new double[data.n];
  data.y = new double[data.n];
  data.log = false;
  data.likelihood = false;

  for (unsigned int i = 0; i < data.n; i++)
  {
    double f = 0.5;
    double a = astro::x_linear
      (i, data.n, a_min+0.05*(1.0-a_min), f+(1.0-f)*a_min);
    data.x[i] = a;
    data.y[i] = a*gsl_pow_2 (a*expansionFunction (a));
  }

  data.model = mf;
  data.dmodel = dm;

  astro::nlModeller nl_model (data);
  double p[2] = {0.3, 0.7};
  nl_model (p);
  q = astro::stdvec (p, 2);
}

double expansionModel::sum_r (double a)
{
  double r = 0.0;
  for (int j = 0; j < n_basis; j++)
    r += c[j]*basis->operator () (a, j);
  return r;
}

double expansionModel::sum_dr (double a)
{
  double dr = 0.0;
  for (int j = 0; j < n_basis; j++)
    dr += gsl_pow_2 (dc[j]*basis->operator () (a, j));
  return sqrt (dr);
}

inline double expansionModel::x_of_a (double a)
{ return (a-a_min)/(1.0-a_min); }

double expansionModel::deltaE (double a)
{
  if (a < a_min)
    return 0.0;
  else
  {
    double x = x_of_a (a);
    double x0 = x_of_a (1.0);
    return fabs(expansionFunction (a)*sum_dr (x)/sum_r (x)-expansionFunction (a)*sum_dr (x0)/sum_r (x0));
  }
}

double expansionModel::expansionFunction_no_norm (double a, void *vp)
{
  if (a < a_min)
    return sqrt ((mf (a, q)+Omega_r0/a)/gsl_pow_3 (a));
  double x = x_of_a (a);
  double r = sum_r (x);
  return delta_a/pow (a_min, exponent)/pow (1.0+x*delta_a, exponent+2.0)/r;
}

double expansionModel::expansionFunction (double a, void *vp)
{
  if (a < a_min)
    return sqrt ((mf (a, q)+Omega_r0/a)/gsl_pow_3 (a));
  double x = x_of_a (a);
  double r = sum_r (x);
  return norm*delta_a/pow (a_min, exponent)/pow (1.0+x*delta_a, exponent+2.0)/r;
}


double expansionModel::dexpansionFunction_da (double a, void * vp)
{
  if (a < a_min)
  {
    double fp = 2.0*q[0]*(a-q[1]);
    return
      0.5*(fp-3.0*mf (a, q)/a-4.0*Omega_r0/a/a)/
      expansionFunction (a)/gsl_pow_3 (a);
  }
  else
  {
    double x = x_of_a (a);
    double r = sum_r (x);
    double rp = 0.0;
    for (int j = 0; j < n_basis; j++)
      rp += c[j]*basis->prime (x, j);
    return
      -expansionFunction (a)/(1.0-a_min)*(rp/r+(exponent+2.0)*delta_a/(1.0+x*delta_a));
  }
}

void expansionModel::print (std::string filename, int n)
{
  std::ofstream output (filename);
  for (int i = 0; i < n; i++)
  {
    double a = astro::x_linear (i, n, a_min, 1.0);
    output
      << a
      << "  " << expansionFunction (a)
      << "  " << deltaE (a)
      << "  " << dexpansionFunction_da (a) << std::endl;
  }
  output.close ();

  std::ofstream output2 ("expansion_basis.d");
  for (int i = 0; i < n; i++)
  {
    double a = astro::x_linear (i, n, a_min, 1.0);
    double x = x_of_a (a);
    output2
      << a;
      for (int j = 0; j < n_basis; j++){
        output2
        << "  " << c[j]*basis->operator () (x, j);}
      output2
      << std::endl;
  }
  output2.close ();
}
