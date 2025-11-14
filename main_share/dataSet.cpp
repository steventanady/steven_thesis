#include <algorithm>
#include <iostream>

#include <astro/utilities/utilities.h>

#include "dataSet.h"

dataSet::dataSet ():
n (0), a_min (1.0), a_max (1.0), delta_a (0.0), normalized (false)
{
  a.clear ();
  d_lum.clear ();
  name.clear ();
  Cov = gsl_matrix_alloc (0, 0);
}

void dataSet::push_back (dataSample& sample)
{
  int n_save = n;
  gsl_matrix * C_save = gsl_matrix_alloc (n_save, n_save);
  gsl_matrix_memcpy (C_save, Cov);

  std::vector<double> a_add = sample.get_a ();
  std::vector<double> d_add = sample.get_distance ();
  std::vector<std::string> name_add = sample.get_name ();
  gsl_matrix * C_add = sample.get_cov_matrix ();

  a.insert (a.end (), a_add.begin (), a_add.end ());
  d_lum.insert (d_lum.end (), d_add.begin (), d_add.end ());
  name.insert (name.end (), name_add.begin (), name_add.end ());

  n = a.size ();

  Cov = gsl_matrix_calloc (n, n);
  for (int j = 0; j < n_save; j++)
    for (int k = 0; k < n_save; k++)
      gsl_matrix_set (Cov, j, k, gsl_matrix_get (C_save, j, k));
  gsl_matrix_free (C_save);

  for (int j = n_save; j < n; j++)
    for (int k = n_save; k < n; k++)
      gsl_matrix_set (Cov, j, k, gsl_matrix_get (C_add, j-n_save, k-n_save));

  a_min = *std::min_element (a.begin (), a.end ());
  a_max = *std::max_element (a.begin (), a.end ());
  delta_a = 0.0;
}

dataSet::dataSet (std::vector<dataSample>& data):
n (0), a_min (1.0), a_max (1.0), delta_a (0.0), normalized (false)
{
  a.clear ();
  d_lum.clear ();
  name.clear ();
  Cov = gsl_matrix_alloc (n, n);

  for (unsigned int i = 0; i < data.size (); i++)
    push_back (data[i]);
}

void dataSet::normalize (double a_min_in)
{
  a_min = a_min_in;
  delta_a = 1.0/a_min-1.0;
  std::vector<double> f (n);
  /**
   * Normalise scale factor and distances.
   */
  for (int i = 0; i < n; i++)
  {
    a[i] = (a[i]-a_min)/(1.0-a_min);
    f[i] = a_min*a_min*(1.0+a[i]*delta_a);
    d_lum[i] *= f[i];
  }
  /**
   * Normalise covariance matrix.
   */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
      double C_ij = gsl_matrix_get (Cov, i, j);
      gsl_matrix_set (Cov, i, j, C_ij*f[i]*f[j]);
    }
  normalized = true;
}

void dataSet::normalize_log ()
{
  std::vector<double> f (n);
  std::vector<double> z (n);
  std::vector<double> scalefactor (n);
  /**
   * Normalise scale factor and distances.
   */
  for (int i = 0; i < n; i++)
  {
    // z[i] = 1./a[i]-1.;
    // double z_min = 1./a_max-1.;
    // double z_max = 1./a_min-1.;
    //
    // a[i] = (log10(z[i])-log10(z_max))/(log10(z_min)-log10(z_max));
    // f[i] = d_lum[i];
    // scalefactor[i]=1./(z_min*(pow(z_min/z_max,a[i]-1.)-1.)+1.);
    // d_lum[i] = log10(scalefactor[i]*d_lum[i]);

    // a[i] = (log10(z[i]+z_min)-log10(z_max+z_min))/(log10(z_min)-log10(z_max+z_min));
    // f[i] = d_lum[i];
    // scalefactor[i]=1./(z_min*(pow(z_min/(z_max+z_min),a[i]-1.)-1.)+1.);
    // d_lum[i] = log10(scalefactor[i]*d_lum[i]);

    a[i] = 1.-log10(a[i])/log10(a_min);
    f[i] = d_lum[i];
    d_lum[i] = pow(a_min,1.-a[i])*d_lum[i];

  }
  /**
   * Normalise covariance matrix.
   */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
      double C_ij = gsl_matrix_get (Cov, i, j);
      gsl_matrix_set (Cov, i, j, C_ij*pow(a_min,1.-a[i])*pow(a_min,1.-a[i]));
    }
  normalized = true;
}


dataSet::~dataSet ()
{ gsl_matrix_free (Cov); }

double dataSet::get_a_min ()
{ return a_min; }

double dataSet::get_a_max ()
{ return a_max; }

double dataSet::get_a (int i)
{ return a[i]; }

double dataSet::get_distance (int i)
{ return d_lum[i]; }

std::vector<double>& dataSet::get_a ()
{ return a; }

std::vector<double>& dataSet::get_distance ()
{ return d_lum; }

std::vector<std::string>& dataSet::get_name ()
{ return name; }

gsl_matrix * dataSet::get_cov_matrix ()
{ return Cov; }

int dataSet::size ()
{ return n; }

bool dataSet::is_normalized ()
{ return normalized; }
