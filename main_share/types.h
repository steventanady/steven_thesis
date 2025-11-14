#ifndef ASTRO_TYPES_H
#define ASTRO_TYPES_H

#include <functional>
#include <vector>

namespace astro
{
  /**
   * \ingroup programming
   *
   * Enumeration for selecting the type of spacing of independent points in
   * an array.
   *
   * \author Felix Fabis, Heidelberg University
   */
  enum spacing_type
  {
    LINEAR_SPACING = 0,
    LOG_SPACING = 1
  };

  /**
   * \ingroup programming
   *
   * Enumeration for the oscillating functions allowed in integrators over
   * oscillating integrals. The structure of the integral kernels is assumed
   * to be \f$f(k,x)w(kx)\f$, where \f$f(k,x)\f$ does not oscillate, but
   * \f$w(kx)\f$ does. The types provided for \f$w(kx)\f$ are
   * \f[\cos(kx)\;,\quad\sin(kx)\;,\quad J_\nu(kx)\;,\quad j_\nu(kx)\;.\f]
   *
   * \author Matthias Bartelmann, modified by Felix Fabis,
   * Heidelberg University , 2016
   */
  enum osc_function_type
  {
    UNSET          = -1,
    OSC_COSINE     =  0,
    OSC_SINE       =  1,
    OSC_TRIG       =  2,
    OSC_CYL_BESSEL =  3,
    OSC_SPH_BESSEL =  4
  };

  /**
   * \ingroup programming
   *
   * Enumeration for the access states to FITS files, allowing the creation
   * of a new FITS file, the read-only access, or the read-write access.
   *
   * \author Matthias Bartelmann, Heidelberg University
   */
  enum fits_access_flag
  {
    FITS_CREATE    = -1,
    FITS_READONLY  = 0,
    FITS_READWRITE = 1
  };

  /**
   * \ingroup programming
   *
   * Enumeration for random-number distributions: uniform or Gaussian
   *
   * \author Matthias Bartelmann, Heidelberg University
   */
  enum random_number_distribution
  {
    RANDOM_UNIFORM = 1,
    RANDOM_GAUSSIAN = 2
  };

  /**
   * \ingroup programming
   *
   * Enumeration for types of CMB spectra: TT, TE, and EE
   *
   * \author Matthias Bartelmann, Heidelberg University
   */
  enum cmb_spectrum_type
  {
    CMB_TT = 1,
    CMB_EE = 2,
    CMB_TE = 3
  };

  /**
   * \ingroup numerics
   *
   * Structure for data sets containing models and model derivatives.
   *
   * \author Matthias Bartelmann, Heidelberg University
   */
  struct data_model
  {
    unsigned int n; // number of data points
    unsigned int p; // number of model parameters
    double * x;     // array of independent data points
    double * y;     // array of dependent data points
    double * dy = NULL; // array of uncertainties in y
    bool log;       // fit model to the log of the data
    bool likelihood;
    std::function<double (double, std::vector<double>)> model;
    std::function<std::vector<double> (double, std::vector<double>)> dmodel;
    std::function<double (double, double, double, std::vector<double>)> mf;
    std::function<std::vector<double> (double, double, double, std::vector<double>)> mJ;
  };

  /**
   * \ingroup numerics
   *
   * Structure containing (r,g,b) values of a pixel.
   *
   * \author Matthias Bartelmann, Heidelberg University
   */
  struct pixel
  { unsigned int red, green, blue; };
}

#endif
