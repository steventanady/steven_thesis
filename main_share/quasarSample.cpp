#include "quasarSample.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <astro/utilities/stringOperations.h>
#include <astro/io/table.h>
#include <astro/io/fits.h>

quasarSample::quasarSample (std::string filename)
{
  std::ifstream test;

  test.open (filename + ".d");
  if (! test.good ())
    throw std::runtime_error ("input file not found");
  test.close ();

  astro::table input_table (filename + ".d");

  std::vector<double> z = input_table.getColumn (0);
  logf_x = input_table.getColumn (1);
  logf_xerr = input_table.getColumn(2);
  logf_uv = input_table.getColumn (3);
  logf_uverr = input_table.getColumn (4);

  n = z.size ();
  a.resize (n);

  for (int i = 0; i < n; i++)
  {
    a[i] = 1.0/(1.0+z[i]);
  }
  a_min = *std::min_element (a.begin (), a.end ());
  a_max = *std::max_element (a.begin (), a.end ());

  for (int i = 0; i<n; i++)
    name.push_back(filename);
}

void quasarSample::normalize(double a_min_in) {

  a_min = a_min_in;

  for (int i = 0; i < n; i++)
  {
    a[i] = (a[i]-a_min)/(1.0-a_min);
  }

}

double quasarSample::get_a_min ()
{ return a_min; }

double quasarSample::get_a_max ()
{ return a_max; }

std::vector<double>& quasarSample::get_a ()
{ return a; }

std::vector<double>& quasarSample::get_logf_x ()
{ return logf_x; }

std::vector<double>& quasarSample::get_logf_xerr ()
{ return logf_xerr; }

std::vector<double>& quasarSample::get_logf_uv ()
{ return logf_uv; }

std::vector<double>& quasarSample::get_logf_uverr ()
{ return logf_uverr; }

std::vector<std::string>& quasarSample::get_name ()
{ return name; }

int quasarSample::size ()
{ return n; }
