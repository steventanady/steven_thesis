#ifndef DATA_SET_H
#define DATA_SET_H

#include <vector>
#include <gsl/gsl_matrix.h>

#include "dataSample.h"

class dataSet
{
protected:
  int n;
  std::vector<double> a, d_lum;
  std::vector<std::string> name;
  gsl_matrix * Cov;
  double a_min, a_max, delta_a;
  bool normalized;

public:
  /**
   * Empty constructor for later addition of data samples.
   */
  dataSet ();

  /**
   * Constructor initialising the data set with a vector of data samples.
   */
  dataSet (std::vector<dataSample>& data);

  /**
   * Destructor, deallocating memory.
   */
  ~dataSet ();

  /**
   * Adds a new data sample to the data set.
   */
  void push_back (dataSample& sample);

  /**
   * Routine normalising the data. This means that that the scale factor is
   * linearly mapped to the interval \f$[0,1]\f$ by applying
   * \f[
   *   a\to x = \frac{a-a_\mathrm{min}}{1-a_\mathrm{min}} =
   *     \frac{a-a_\mathrm{min}}{a_\mathrm{min}\delta a}
   * \f]
   * to the scale factor and
   * \f[
   *   D_\mathrm{L}\to d_\mathrm{L} =
   *     a_\mathrm{min}^2(1+x\delta a)D_\mathrm{L}
   * \f]
   * with
   * \f[
   *   \delta = \frac{1-a_\mathrm{min}}{a_\mathrm{min}}\;.
   * \f]
   * The elements of the covariance matrix are scaled by
   * \f[
   *   C_{ij}\to C_{ij}f_if_j
   * \f]
   * with
   * \f[
   *   f_i = a_\mathrm{min}^2(1+x_i\delta a)\;.
   * \f]
   */
  void normalize (double a_min_in);
  void normalize_log ();

  double get_a_min ();
  double get_a_max ();

  double get_a (int i);
  double get_distance (int i);

  std::vector<double>& get_a ();
  std::vector<double>& get_distance ();
  std::vector<std::string>& get_name ();

  gsl_matrix * get_cov_matrix ();

  int size ();

  bool is_normalized ();
};

#endif
