#include <string>
#include <iterator>
#include <cmath>

#include <astro/utilities/shiftedChebyshev.h>
#include <astro/utilities/stringOperations.h>
#include <astro/utilities/integrator.h>
#include <astro/utilities/deSolver.h>
#include <astro/io/clArguments.h>
#include <astro/io/functionWriter.h>

#include "dataSet.h"
#include "quasarSample.h"
#include "distanceModel.h"
#include "distanceModes.h"
#include "expansionModel.h"
#include "friedmannFit.h"
#include "growthModels.h"


/**
 * Main routine to demonstrate the operation of the code pieces. Reads
 * type-Ia supernova (and/or bao, qso, cmb) data from file, constructs a
 * distance model from it, converts this to a model for the expansion
 * function and finally uses that to calculate the growth factor.
 */


int main (int argc, char * argv[])
{
  /**
   * Read command-line arguments
   */
  astro::clArguments ca (argc, argv, "help.txt");
  if (ca.set ("-h"))
  {
    ca.print_help ();
    return 1;
  }

  std::string input = ca.get ("-i", "none");
  std::string qinput = ca.get ("-q", "none");

  int n_basis = ca.get ("-b", 4);

  /**
   * Set exponent in the expression E (a) = a**(n)*R (a),
   * the function R (a) is then modelled with polynomials
   */
  double exponent = ca.get ("-e", 0.0);

  /**
   * Read input file names and offsets from argument list
   */
  std::vector<std::string> list = astro::split (input, ";");

  /**
   * Read data sample(s) from file(s) and normalize the scale factor to
   * the interval [0,1] if needed
   */

  dataSet * sample;

  if (input != "none"){

    sample = new dataSet;

    for (unsigned int i = 0; i < list.size (); i++)
    {
      dataSample dsample ("data/"+list[i]);
      sample->push_back (dsample);
    }
  }
  else
    sample = NULL;


  double a_min;
  quasarSample * qsample;

  if (qinput!="none"){
    qsample = new quasarSample ("data/"+qinput);
    if (input != "none")
      a_min = std::min (sample->get_a_min (),qsample->get_a_min ());
    else
      a_min = qsample->get_a_min ();
    qsample->normalize (a_min);
  }
  else {
    a_min = sample->get_a_min();
    qsample = NULL;
  }

  if (input != "none")
    sample->normalize (a_min);


  //const double Omega_m0 = 0.3;

  /**
   * Initialize the shifted Chebyshev polynomials on the interval [0,1],
   * write them to a file and initialize the distance and expansion
   * modes
   */
  astro::shiftedChebyshev basis;
  basis.setParameters (n_basis, 0.0, 1.0);
  basis.print ("basis.d");

  distanceModes d_modes (&basis, a_min, exponent);
  d_modes.setParameters (n_basis, 0.0, 1.0);

  /**
   * Write distance and expansion-function modes to files
   */
  d_modes.print ("d_modes.d");

  /**
   * Construct the distance model by adapting the basis to the
   * sample
   */

  distanceModel distance (sample, qsample, &d_modes, exponent);

  distance.print ("distance.d");
  distance.report ("coefficients.d", input+";"+qinput, a_min, exponent);

  std::vector<double> c = distance.get_c ();
  std::vector<double> dc = distance.get_delta_c ();

  std::cout<<"Finished calculating Distance Model"<<std::endl;

  /**
   * Calculates expansion model from the distance model
   */
  expansionModel expansion (&basis, c, dc, a_min, exponent);
  expansion.print ("expFunction.d");

  std::cout<<"Finished calculating Expansion Function"<<std::endl;

  /**
   * Fits spatially-flat Friedmann model to the expansion function
   */
  friedmannFit fit ("expFunction.d", &basis, c.size ());
  fit.report ("coefficients.d");
  fit.de_model ("darkEnergy.d");

  double Omega_m0 = fit.get_omega_red ();

  std::cout<<"Finished Friedmann Fit"<<std::endl;

  /**
   * Calculates the growth factor from the expansion
   * function and its derivative
   */
  growthModels growth (&expansion, a_min, 0.35);
  growth.print ("growthFunction.d");

  std::cout<<"Finished calculating growth model"<<std::endl;

  return 0;
}
