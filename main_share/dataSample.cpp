#include "dataSample.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <astro/utilities/stringOperations.h>
#include <astro/io/table.h>
#include <astro/io/fits.h>

dataSample::dataSample (std::string filename)
{
  std::ifstream test;

  test.open (filename + ".d");
  if (! test.good ())
    throw std::runtime_error ("input file not found");
  test.close ();

  astro::table input_table (filename + ".d");
  double offset = input_table.read_key ("offset");
  double so = input_table.read_key ("delta_offset");

  std::vector<double> z = input_table.getColumn (0);
  mu = input_table.getColumn (1);
  std::vector<double> sp = input_table.getColumn (2);
  std::vector<double> sm = input_table.getColumn (3);
  n = z.size ();
  a.resize (n);
  d_lum.resize (n);
  double cf = 0.2*M_LN10;
  for (int i = 0; i < n; i++)
  {
    a[i] = 1.0/(1.0+z[i]);
    d_lum[i] = exp (M_LN10*(0.2*(offset+mu[i])+1.0))/1.0e6/3.0e3;
  }
  a_min = *std::min_element (a.begin (), a.end ());
  a_max = *std::max_element (a.begin (), a.end ());

  for (int i = 0; i<n; i++)
    name.push_back(filename);

  Cov = gsl_matrix_alloc (n, n);

  test.open (filename + "_CM.fits");
  if (! test.good ())
  {
    gsl_matrix_set_zero (Cov);
    for (int i = 0; i < n; i++)
    {
      double sigsqr =
        gsl_pow_2 (cf*d_lum[i])*(gsl_pow_2 (sp[i]+sm[i])+gsl_pow_2 (so));
      gsl_matrix_set (Cov, i, i, sigsqr);
    }
  }
  else
  {
    astro::fits inputCM (filename + "_CM.fits", astro::FITS_READONLY);
    std::vector<long> dim = inputCM.get_dim ();
    if (dim[0] != n)
      std::cerr
        << "data size mismatch, ignoring covariance matrix" << std::endl;
    std::vector<double> data;
    inputCM.read (data);

    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
      {
        double C_ij = cf*cf*data[i*n+j]*d_lum[i]*d_lum[j];
        gsl_matrix_set (Cov, i, j, C_ij);
      }
  }
  test.close ();
}

double dataSample::get_a_min ()
{ return a_min; }

double dataSample::get_a_max ()
{ return a_max; }

double dataSample::get_a (int i)
{ return a[i]; }

double dataSample::get_distance (int i)
{ return d_lum[i]; }

std::vector<double>& dataSample::get_a ()
{ return a; }

std::vector<double>& dataSample::get_distance ()
{ return d_lum; }

std::vector<std::string>& dataSample::get_name ()
{ return name; }

gsl_matrix * dataSample::get_cov_matrix ()
{ return Cov; }

int dataSample::size ()
{ return n; }
