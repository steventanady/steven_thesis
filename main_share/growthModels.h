#ifndef GROWTH_MODELS_H
#define GROWTH_MODELS_H

#include <astro/utilities/integrator.h>
#include <astro/utilities/functionTable.h>
#include <astro/cosmology/abstractGrowthFactor.h>

#include "expansionModel.h"

class growthModels: public astro::abstractGrowthFactor
{
protected:
  expansionModel * exp_model;
  double a_min;
  double omega_m0;

  astro::functionTable D_of_a;
  astro::functionTable a_of_D;
  astro::functionTable growth_tab;
  astro::functionTable deltaDplus_tab; 
public:
  /**
   * Constructor initializing the growth model with an expansion
   * model together with the lower limit of the scale factor, matter
   * density parameter today and the exponent \c n in a local approximation
   * of the growth factor \f$D_+\f$ by a power law, \f$D_+ = a^n\f$.
   */
  growthModels
    (expansionModel * exp_model_in, double a_min_in,
     double omega_m0_in);

  //~growthModels ();

 /**
  * Returns the growth factor over the scale factor \c a at a given \c a
  */
  double operator () (double a, double k = 0.0);

 /**
  * Returns the growth factor at a given scale factor \c a
  */
  double Dplus (double a);

  double Dplus_marg (double a, double mean, double sigma);

  double Dplus_marg_std (double a, double mean, double sigma);

  double deltaDplus_marg (double a, double mean, double sigma);

  double gaussFilter (double x, double mu, double sigma);

  double dlnDplus_dlna (double a);

  double dDplus_dOmega_m0 (double a);

  double dDplus_dOmega_m02 (double a);

  double index (double a, double Omega_m0);

  double growth (double a,double Omega_m0);

 /**
  * Returns the reduced growth factor at a given scale factor \c a
  * using the differential equation for the reduced growth
  * \f$G_+ (a) = D_+ (a)/a\f$
  */
  double reducedGrowth (double a, double Omega_m0);

 /**
  * Returns a factor to multiply the gravitational constant with.
  */
  double gravitationalConstant (double a, double k = 1.0);

  /**
   * Derivative of the growth factor \c D_+ with respect to \c E
   */
  double dDplus_dE (double x, double a, double Omega_m0);

  /**
   * Returns the uncertainty on the expansion function
   */
  double deltaE (double a);

  /**
   * Returns the uncertainty of the growth factor
   */
  double deltaDplus (double a, double Omega_m0);

 /**
  * Writes the growth factor to a file with the given name,
  * evaluating it at the number of points given.
  */
  void print (std::string filename, int n = 128);
};

#endif
