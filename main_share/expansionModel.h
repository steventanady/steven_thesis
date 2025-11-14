#ifndef EXP_MODEL_H
#define EXP_MODEL_H

#include <astro/utilities/functionSet.h>
#include <astro/cosmology/expansionModel.h>

class expansionModel: public astro::expansionModel
{
protected:
  astro::functionSet * basis;
  int n_basis;
  std::vector<double> c, dc, q;
  double a_min, delta_a;
  double exponent;
  double norm;
  double Omega_r0;

  /**
   * Returns the inverse expansion function \f$e = (aE)^{-1}\f$ as a sum over
   * modes
   */
  double sum_r (double a);

  /**
   * Returns the uncertainty on the inverse expansion function as a sum over
   * modes and uncertainties in distance coefficients
   */
  double sum_dr (double a);

  /**
   * Returns the normalised scale factor \c x at the scale factor \c a.
   */
  inline double x_of_a (double a);

  static double mf (double a, std::vector<double> p);
  static std::vector<double> dm (double a, std::vector<double> p);
  void fitModel ();
  void init_model ();
public:
  /**
   * Constructor initializing the expansion model with a distance model
   * together with the lower limit of the scale factor and exponent
   * of the prefactor.
   */
  expansionModel
      (astro::functionSet * basis_in, std::vector<double>& c_in,
       std::vector<double>& dc_in, double a_min_in, double exponent_in);


  /**
   * Returns the expansion function at a given scale factor \c a
   */
  double expansionFunction_no_norm (double a, void * vp = NULL);



  /**
   * Returns the expansion function at a given scale factor \c a
   */
  double expansionFunction (double a, void * vp = NULL);

  /**
   * Returns the uncertainty on the expansion function as a function of the
   * scale factor \c a.
   */
  double deltaE (double a);

  /**
   * Returns the derivative of the expansion function with respect to
   * the scale factor at a given scale factor \c a
   */
  double dexpansionFunction_da (double a, void * vp = NULL);

  /**
   * Write the expansion function and its derivative to a file with the
   * given name, evaluating it at the number of points given.
   */
  void print (std::string filename, int n = 1028);

};

#endif
