#ifndef FRIEDMANN_FIT_H
#define FRIEDMANN_FIT_H

#include <string>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <astro/utilities/types.h>
#include <astro/utilities/functionSet.h>

class friedmannFit
{
protected:
  double a_min, delta_a;
  double Omega_m0, delta_Omega_m0;
  double Omega_m0_red, delta_Omega_m0_red;
  astro::data_model sample, sample_red;
  unsigned int n_fit;
  std::vector<double> x, y, d;
  std::vector<double> c_fit;

  double e (double a);

  static double E (double a, std::vector<double> p);

  static std::vector<double> dE (double a, std::vector<double> p);

  void init_coeffs (astro::functionSet * basis);
public:
  friedmannFit
    (std::string filename, astro::functionSet * basis, unsigned int n_fit_in);

  double significance (std::vector<double>& c_mean, gsl_matrix * F);

  void de_model (std::string filename);

  void report (std::string filename);

  double get_omega ();

  double get_omega_red ();
};

#endif
