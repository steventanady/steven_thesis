#ifndef DISTANCE_MODEL_H
#define DISTANCE_MODEL_H

#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <astro/utilities/utilities.h>
#include <astro/utilities/functionSet.h>
#include "dataSet.h"
#include "quasarSample.h"

/**
 * Creates a model for the luminosity distance from a measured sample
 * and a function set used as a basis.
 */
class distanceModel
{
protected:
  int n_data, n_basis, n_qso, n_dataq;
  double a_min, a_max, a_min_s, exponent;
  bool nsample, quasars;
  /**
   * Eigenvalues of the Fisher matrix, the Fisher matrix itself, and the
   * rotation matrix diagonalizing it:
   */
  gsl_vector * v;
  gsl_matrix * R;
  /**
   * Auxiliary matrices needed during the calculations:
   */
  gsl_matrix * P, * Cov, *tCov, * invtCov, * J, * JTinvCov, * F, * cov, * invcov;
  /**
   * Vectors for the coefficients and their uncertainties:
   */
  std::vector<double> p, dp, c, dc;

  std::vector<double> scale, distance, scaleq, logf_x, logf_xerr, logf_uv,
  logf_uverr;
  std::vector<std::string> name, nameq;
  astro::functionSet * basis;

  struct data_model
  {
    int n; // number of data points
    int nd;
    int nq;
    int p; // number of model parameters
    int q;
    double * x;     // array of independent data points
    double * y;     // array of dependent data points
    double * dy = NULL; // array of uncertainties in y
    double * z = NULL; //second array of dependent data points
    double * dz = NULL; //array of uncertainties in z
    gsl_matrix * cov = NULL; //covariance matrix
    gsl_matrix * invcov = NULL; //covariance matrix
    std::string * type = NULL; // array of source type
    double a_min = 0.0;
    astro::functionSet * basis = NULL;
  };

  double mf (std::string s, double x, double y, double dy, double z, double dz,
     std::vector<double> p);

  std::vector<double> mJ (std::string s, double x, double y,
     double dy, double z, double dz, std::vector<double> p);

  static double my_f (const gsl_vector *v, void *d);

  static void my_df (const gsl_vector *v, void *params, gsl_vector *df);

  static void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df);

  /**
   * Calculate Fisher-Matrix
   */
  void init_Fisher (data_model dm);

  /**
   * Calculates the basis-function coefficients
   */
  void init_Coeffs ();

  /**
   * Check whether all coefficients are significant and set insignificant
   * coefficients to zero
   */
  bool clean_Coeffs ();
public:
  /**
   * Constructor initializing the distance model with a measurement sample
   * and with a pointer to a function set
   */
  distanceModel (dataSet * sample, quasarSample * qsample,
    astro::functionSet * basis_in, double exponent_in);

  /**
   * Destructor releasing memory allocated for GSL matrices
   */
  ~distanceModel ();

  /**
   * Returns the model distance to the scale factor \c a
   */
  double operator () (double a);

  double error (double a);

  /**
   * Returns the derivative of the distance to the scale factor \c a
   */
  double prime (double a);

  /**
   * Returns the second derivative of the distance to the scale factor \c a
   */
  double pprime (double a);

  /**
   * Returns the contribution \f$e_i\f$ of the basis functions \f$p_i\f$ to
   * the expansion function.
   */
  double e (double a, int i);

  /**
   * Returns the contribution \f$e_i'\f$ of the derivative \f$p_i'\f$ of the
   * basis function \f$p_i\f$.
   */
  double e_prime (double a, int i);

  /**
   * Returns the coefficient \f$c_i\f$ of the \f$i\f$-th basis function
   */
  double get_c (int i);

  /**
   * Returns the vector of all basis-function coefficients
   */
  std::vector<double> get_c ();

  /**
   * Returns the uncertainty on the \f$i\f$-th coefficient
   */
  double get_delta_c (int i);

  /**
   * Returns the vector of all uncertainties
   */
  std::vector<double> get_delta_c ();

  /**
   * Returns the Fisher matrix for the coefficients
   */
  gsl_matrix * get_covariance_matrix ();

  /**
   * Write the distance model to file, together with the measured distances.
   */
  void print (std::string filename, int n = 128);

  /**
   * Write the expansion coefficients, the Fisher matrix and other
   * information to the file with the given name.
   * \param filename name of the output file
   * \param input string characterizing the input data
   * \param a_min minimum scale factor in data sample
   */
  void report (std::string filename, std::string input, double a_min, double exponent);
};

#endif
