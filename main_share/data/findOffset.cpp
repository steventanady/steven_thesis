/**
 * g++ -O -I$HOME/.local/include -o findOffset.x findOffset.cpp -L$HOME/.local/lib -lastrocpp -lgsl -lgslcblas
 */
#include <astro/io/table.h>
#include <astro/utilities/utilities.h>
#include <iostream>
#include <algorithm>
#include <gsl/gsl_statistics_double.h>

/**
 * Structure for a data set containing an independent variable \c x, an
 * dependent variable \c y, a standard deviation \c s, and the number \c n
 * of data points.
 */
struct dataset
{
  int n;
  std::vector<double> x, y, s;
};

/**
 * Routine returning two vectors containing the mean and the variance of a
 * data set in given bins.
 */
std::vector<std::vector<double> > findMean
  (dataset& data, std::vector<double>& bins)
{
  int n_bins = bins.size ()-1;
  std::vector<std::vector<double> > mv (2);
  for (int i = 0; i < n_bins; i++)
  {
    std::vector<double> x, w;
    for (int j = 0; j < data.n; j++)
    {
      if (data.x[j] >= bins[i] && data.x[j] < bins[i+1])
      {
        x.push_back (data.y[j]);
        w.push_back (1.0/gsl_pow_2 (data.s[j]));
      }
    }
    mv[0].push_back (gsl_stats_wmean (w.data (), 1, x.data (), 1, x.size ()));
    mv[1].push_back
      (gsl_stats_wvariance (w.data (), 1, x.data (), 1, x.size ()));
  }
  return mv;
}

int main ()
{
  /**
   * Open two data files and read data sets from them.
   */
  dataset d_1, d_2;
  astro::table tab_1 ("pantheonSample.d");
  astro::table tab_2 ("baoSample.d");

  double offset = tab_1.read_key ("offset");

  d_1.x = tab_1.getColumn (0);
  d_1.y = tab_1.getColumn (1);
  d_1.s = tab_1.getColumn (2);
  d_1.n = d_1.x.size ();

  d_2.x = tab_2.getColumn (0);
  d_2.y = tab_2.getColumn (1);
  d_2.s = tab_2.getColumn (2);
  d_2.n = d_2.x.size ();

  /**
   * Prepare bins such that each bin contains the same number \c delta_n
   * of data points.
   */
  int delta_n = 16;
  std::vector<double> s = d_1.x;
  std::sort (s.begin (), s.end ());
  std::vector<double> bins;
  for (int i = 0; i < s.size (); i += delta_n)
    bins.push_back (s[i]);

  /**
   * Determine means and variances of the two data sets in these bins.
   */
  std::vector<std::vector<double> > mv_1 = findMean (d_1, bins);
  std::vector<std::vector<double> > mv_2 = findMean (d_2, bins);
  std::vector<std::vector<double> > mv_3 (2);

  /**
   * In those bins containing means and variances for both data sets,
   * determine the offset and their statistical weights.
   */
  for (int i = 0; i < mv_2[1].size (); i++)
    if (mv_2[0][i] > 0 && mv_2[1][i] > 0)
    {
      mv_3[0].push_back (mv_2[0][i]-mv_1[0][i]);
      mv_3[1].push_back (1.0/(mv_2[1][i]+mv_1[1][i]));
    }

  /**
   * Finally, determine and print the offset averaged over all bins and its
   * uncertainty.
   */
  std::cout
    << "# offset = "
    << offset-gsl_stats_wmean
      (mv_3[1].data (), 1, mv_3[0].data (), 1, mv_3[0].size ())
    << std::endl << "# delta_offset = "
    << sqrt (gsl_stats_wvariance
      (mv_3[1].data (), 1, mv_3[0].data (), 1, mv_3[0].size ()))
    << std::endl << "#" << std::endl;

  return 0;
}
