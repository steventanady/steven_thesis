// g++ -O -I$HOME/.local/include -o covMatrix covMatrix.cpp -L$HOME/.local/lib -lastrocpp -lgsl -lgslcblas

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>

#include <astro/io/fits.h>
#include <astro/io/clArguments.h>
#include <astro/utilities/stringOperations.h>

int main (int argc, char * argv[])
{
  astro::clArguments ca (argc, argv);
  std::string input_file = ca.get ("-i", "pantheon+Sample_CM");

  std::string line;
  std::vector<long> dim (2, 0);
  std::vector<double> data;
  std::ifstream input (input_file + ".d");
  if (!input.good ())
    throw std::runtime_error ("input file not found");
  while (std::getline (input, line))
  {
    if (! astro::is_empty (line))
    {
      std::istringstream buffer (line);
      std::vector<double> row
        ((std::istream_iterator<double> (buffer)),
        std::istream_iterator<double> ());
      if (dim[0] == 0)
        dim[0] = row.size ();
      dim[1]++;
      data.insert (data.end (), row.begin (), row.end ());
    }
  }
  input.close ();

  if (dim[0] != dim[1])
    throw std::length_error ("numbers of rows and columns differ");
  if (int (data.size ()) != dim[0]*dim[1])
    throw std::length_error ("data size mismatch");

  astro::fits output (input_file + ".fits", astro::FITS_CREATE);
  output.push_back (data, dim);
  output.close ();

  return 0;
}
