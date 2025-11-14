#include "friedmannFit.h"

#include <iostream>
#include <algorithm>
#include <astro/io/table.h>
#include <astro/utilities/utilities.h>
#include <astro/utilities/nlModeller.h>

friedmannFit::friedmannFit
  (std::string filename, astro::functionSet * basis, unsigned int n_fit_in):
Omega_m0 (0.5), n_fit (n_fit_in)
{
  astro::table data (filename);
  x = data.getColumn (0);
  y = data.getColumn (1);
  d = data.getColumn (2);

  a_min = *std::min_element (x.begin (), x.end ());
  delta_a = 1.0/a_min-1.0;

  sample.n = x.size ()-1;
  sample.p = 1;
  sample.x = new double[sample.n];
  sample.y = new double[sample.n];
  sample.dy = new double[sample.n];
  for (unsigned int i = 0; i < (sample.n); i++)
  {
    sample.x[i] = x[i];
    sample.y[i] = y[i];
    sample.dy[i] = d[i];
  }
  sample.log = false;
  sample.model  = E;
  sample.dmodel = dE;

  astro::nlModeller model (sample);
  double p[1] = {Omega_m0};
  model (p);
  gsl_matrix * cov = gsl_matrix_alloc (sample.p, sample.p);
  model.get_cov_matrix (p, cov);
  Omega_m0 = p[0];
  delta_Omega_m0 = sqrt (gsl_matrix_get (cov, 0, 0));

  init_coeffs (basis);

  //fit Omega_m0 for small redshifts only

  std::vector<double> x_red;
  std::vector<double> y_red;
  std::vector<double> d_red;

  for(int i=0; i < x.size(); i++){

    if (x[i] >= 0.4 & x[i] != 1.){

      x_red.push_back(x[i]);
      y_red.push_back(y[i]);
      d_red.push_back(d[i]);
    }
  }

  sample_red.n = x_red.size ();
  sample_red.p = 1;
  sample_red.x = new double[sample_red.n];
  sample_red.y = new double[sample_red.n];
  sample_red.dy = new double[sample_red.n];
  for (unsigned int i = 0; i < (sample_red.n); i++)
  {
    sample_red.x[i] = x_red[i];
    sample_red.y[i] = y_red[i];
    sample_red.dy[i] = d_red[i];
  }
  sample_red.log = false;
  sample_red.model  = E;
  sample_red.dmodel = dE;

  astro::nlModeller model_red (sample_red);
  double p_red[1] = {0.3};
  model_red (p_red);
  gsl_matrix * cov_red = gsl_matrix_alloc (sample_red.p, sample_red.p);
  model_red.get_cov_matrix (p_red, cov_red);
  Omega_m0_red = p_red[0];
  delta_Omega_m0_red = sqrt (gsl_matrix_get (cov_red, 0, 0));

}

double friedmannFit::e (double x)
{
  std::vector<double> p (1, Omega_m0);
  double a = a_min*(1.0+delta_a*x);
  return delta_a/E (a, p)/gsl_pow_2 (1.0+delta_a*x);
}

double friedmannFit::E (double a, std::vector<double> p)
{ return sqrt (p[0]/a/a/a+(1.0-p[0])); }

std::vector<double> friedmannFit::dE (double a, std::vector<double> p)
{
  std::vector<double> v (1, 0.5/E (a, p)*(1.0/a/a/a-1.0));
  return v;
}

void friedmannFit::init_coeffs (astro::functionSet * basis)
{
  c_fit.resize (n_fit);
  for (unsigned int i = 0; i < n_fit; i++)
    c_fit[i] = basis->scalar_product
      ([&] (double x) { return e (x); }, i);
}

double friedmannFit::significance (std::vector<double>& c_mean, gsl_matrix * F)
{
  if (n_fit == 0)
    throw std::runtime_error ("coefficients undetermined");
  if (c_mean.size () != n_fit)
    throw std::length_error ("coefficient vectors are of different size");

  gsl_vector * dc = gsl_vector_alloc (n_fit);
  for (unsigned int i = 0; i < n_fit; i++)
    gsl_vector_set (dc, i, c_fit[i]-c_mean[i]);
  gsl_vector * Fdc = gsl_vector_alloc (n_fit);
  gsl_blas_dgemv (CblasNoTrans, 1.0, F, dc, 1.0, Fdc);
  double r;
  gsl_blas_ddot (dc, Fdc, &r);
  return sqrt (r/double (n_fit));
}

void friedmannFit::de_model (std::string filename)
{
  double Omega_q0 = 1.0-Omega_m0;
  std::vector<double> p (1, Omega_m0);
  std::ofstream output (filename);
  for (unsigned int i = 0; i < x.size (); i++)
  {
    double e_1 = y[i];
    double e_2 = E (x[i], p);
    output
      << x[i] << "  "
      << (e_1*e_1-e_2*e_2)/Omega_q0+1.0 << "  "
      << 2.0*e_2*d[i]/Omega_q0
      << std::endl;
  }
  output.close ();
}

void friedmannFit::report (std::string filename)
{
  std::ofstream output (filename, std::ios_base::app);
  output
    << std::endl
    << "Results of Friedmann fit to expansion function:" << std::endl
    << "Omega_m0 fit:   " << Omega_m0 << std::endl
    << "delta Omega_m0: " << delta_Omega_m0 << std::endl
    << "Omega_m0_red fit:   " << Omega_m0_red << std::endl
    << "delta Omega_m0_red: " << delta_Omega_m0_red << std::endl
    << "fit coefficients: "
    << astro::str (c_fit) << std::endl;
  output.close ();
}

double friedmannFit::get_omega ()
{ return Omega_m0; }

double friedmannFit::get_omega_red ()
{ return Omega_m0_red; }
