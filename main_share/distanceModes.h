#ifndef DISTANCE_MODES_H
#define DISTANCE_MODES_H

#include <astro/utilities/functionSet.h>

class distanceModes: public astro::functionSet
{
protected:
  astro::functionSet * basis;
  
  double a_min;
  double exponent;
  double delta_a;
  
public:
  distanceModes
  (astro::functionSet * basis_in, double a_min_in, double exponent_in);
  
  double operator () (double a, int n);
  
  double prime (double a, int n);
  
  double pprime (double a, int n);
};

#endif
