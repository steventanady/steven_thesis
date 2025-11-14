#ifndef QUASAR_SAMPLE_H
#define QUASAR_SAMPLE_H

#include <vector>
#include <string>

#include <gsl/gsl_matrix.h>

/**
 * \class dataSample
 *
 * \brief Class providing access to quasar fluxes stored in a file.
 *
 * Class defining methods to access a table of distance-modulus measurements
 * stored in a file. The class assumes that the columns in the table are
 * arranged like this:
 *   column 1: redshift
 *   column 2: log x-ray flux
 *   column 3: deviation of log x-ray flux
 *   column 4: log uv flux
 *   column 5: deviation of log uv flux
 *
 * \author Alexander Oestreicher, Heidelberg University
 */
class quasarSample
{
  int n; // number of data points found in the sample
  double a_min, a_max; // minimum and maximum scale factors

  std::vector<double>
    a,     // scale factors
    logf_x, logf_xerr, logf_uv, logf_uverr;

  std::vector<std::string> name; //sample name

public:
  quasarSample (std::string filename);

  void normalize (double a_min_in);

  double get_a_min ();
  double get_a_max ();

  std::vector<double>& get_a ();
  std::vector<double>& get_logf_x ();
  std::vector<double>& get_logf_xerr ();
  std::vector<double>& get_logf_uv ();
  std::vector<double>& get_logf_uverr ();
  std::vector<std::string>& get_name ();

  int size ();
};

#endif
