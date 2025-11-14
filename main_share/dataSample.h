#ifndef DATA_SAMPLE_H
#define DATA_SAMPLE_H

#include <vector>
#include <string>

#include <gsl/gsl_matrix.h>

/**
 * \class dataSample
 *
 * \brief Class providing access to distance-modulus data stored in a file.
 *
 * Class defining methods to access a table of distance-modulus measurements
 * stored in a file. The class assumes that the columns in the table are
 * arranged like this:
 *   column 1: redshift
 *   column 2: distance modulus
 *   column 3: positive deviation of distance modulus
 *   column 4: negative deviation of distance modulus
 *
 * \author Matthias Bartelmann, Sophia Haude, Heidelberg University
 */
class dataSample
{
  int n; // number of data points found in the sample
  double a_min, a_max; // minimum and maximum scale factors
//   double delta_a;      // used for normalising the data, 0 by default
//   bool normalize;
  std::vector<double>
    a,     // scale factors
    mu,    // distance moduli
    d_lum, // luminosity distances in c/H_0
    sigma; // uncertainty of luminosity distances in c/H_0
  std::vector<std::string> name; //sample name
  gsl_matrix * Cov; // pointer to covariance matrix

public:
  dataSample (std::string filename);

  double get_a_min ();
  double get_a_max ();

  double get_a (int i);
  double get_distance (int i);

  std::vector<double>& get_a ();
  std::vector<double>& get_distance ();
  std::vector<std::string>& get_name ();

  gsl_matrix * get_cov_matrix ();

  int size ();
};

#endif
