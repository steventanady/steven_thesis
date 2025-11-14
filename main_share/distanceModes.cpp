#include "distanceModes.h"

#include <iostream>
#include <astro/utilities/integrator.h>

distanceModes::distanceModes
  (astro::functionSet * basis_in, double a_min_in, double exponent_in):
basis (basis_in), a_min (a_min_in), exponent (exponent_in)
{
  n_basis = basis->size ();
  delta_a = 1.0/a_min - 1.0;
}

double distanceModes::operator () (double a, int n)
{
  astro::integrator integrate
    ([&] (double x)
    { return basis->phi (x, n)*pow (1.0+x*delta_a, exponent); },
    1.0e-4, 1.0e-4);
  return integrate (a, 1.0);
}

double distanceModes::prime (double a, int n)
{ return -basis->phi (a, n)*pow (1.0+a*delta_a, exponent); }

double distanceModes::pprime (double a, int n)
{
  return prime (a, n) * (basis->phi_prime (a, n)/basis->phi (a, n)
         + exponent*delta_a/(1.0+a*delta_a));
}
