#ifndef ASTRO_NL_MODELLER_H
#define ASTRO_NL_MODELLER_H

// #include <gsl/gsl_version.h>

#if GSL_MAJOR_VERSION > 1
#include <gsl/gsl_multifit_nlinear.h>
#else
#include <gsl/gsl_multifit_nlin.h>
#endif
#include "types.h"

namespace astro
{
  /**
   * \ingroup numerics
   * 
   * \class nlModeller
   * 
   * \brief Provides methods for adapting non-linear models to a data set
   * 
   * Class providing methods for adapting non-linear models to a data set,
   * based on GSL functions.
   * 
   * \author Matthias Bartelmann, Heidelberg University
   */
  class nlModeller
  {
  protected:
    astro::data_model data;
#if GSL_MAJOR_VERSION > 1
    const gsl_multifit_nlinear_type * type = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_parameters parameters =
      gsl_multifit_nlinear_default_parameters ();
    gsl_multifit_nlinear_workspace * work;
    gsl_multifit_nlinear_fdf fdf;
#else
    const gsl_multifit_fdfsolver_type * type = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_function_fdf fdf;
#endif    
  
    static int mf (const gsl_vector * x, void * d, gsl_vector * f);
    static int mJ (const gsl_vector * x, void * d, gsl_matrix * J);
#if GSL_MAJOR_VERSION > 1
#else
    static int fJ (const gsl_vector * x, void * d, gsl_vector * f, gsl_matrix * J);
#endif
  public:
    /**
     * Initialize the non-linear modeller with a data model
     * 
     * \param data_in input data model
     * \param jacobian flag indicating whether Jacobian should be used in fit
     */
    nlModeller (astro::data_model& data_in, bool jacobian = true);
  
    /**
     * Returns the values of the best-fitting parameters
     * 
     * \param q array of the best-fitting parameters
     * \param n_max maximum number of iterations
     * \param epsabs absolute error bound
     * \param epsrel relative error bound
     */
    void operator ()
      (double * q, const unsigned int n_max = 100,
       const double epsabs = 1.0e-4, const double epsrel = 1.0e-4);
    
    /**
     * Returns the covariance matrix of the best-fitting parameters
     * 
     * \param p array of the best-fitting parameters
     * \param cov pointer to the covariance matrix
     * \param epsrel relative error bound
     */
    void get_cov_matrix
      (double * p, gsl_matrix * cov, const double epsrel = 1.0e-4);
  };
}
#endif
