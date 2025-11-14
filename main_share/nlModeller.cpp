#include "nlModeller.h"

#include <vector>
#include "utilities.h"

int astro::nlModeller::mf (const gsl_vector * x, void * d, gsl_vector * f)
{
  astro::data_model dm = *reinterpret_cast<astro::data_model *> (d);
  std::vector<double> vx = astro::stdvec (x);
  if (dm.likelihood)
  {
    for (unsigned int i = 0; i < dm.n; i++)
    {
      gsl_vector_set (f, i, dm.mf (dm.x[i], dm.y[i], dm.dy[i], vx));
    }
  }
  else
  {
    for (unsigned int i = 0; i < dm.n; i++)
    {
      double y = (dm.log) ? log (dm.y[i]) : dm.y[i];
      double dy;
      if (dm.dy)
        dy = (dm.log) ? fabs (dm.dy[i]/dm.y[i]) : fabs (dm.dy[i]);
      else
        dy = 1.0;
      gsl_vector_set (f, i, (dm.model (dm.x[i], vx)-y)/dy);
    }
  }
  return GSL_SUCCESS;
}

int astro::nlModeller::mJ (const gsl_vector * x, void * d, gsl_matrix * J)
{
  astro::data_model dm = *reinterpret_cast<astro::data_model *> (d);
  std::vector<double> vx = astro::stdvec (x);
  if (dm.likelihood)
  {
    for (unsigned int i = 0; i < dm.n; i++)
    {
      std::vector<double> mJ = dm.mJ (dm.x[i], dm.y[i], dm.dy[i], vx);
      for (unsigned int j = 0; j < df.size (); j++)
        gsl_matrix_set (J, i, j, mJ[j]);
    }
  }
  else
  {
    for (unsigned int i = 0; i < dm.n; i++)
    {
      std::vector<double> df = dm.dmodel (dm.x[i], vx);
      double dy;
      if (dm.dy)
        dy = (dm.log) ? fabs (dm.dy[i]/dm.y[i]) : fabs (dm.dy[i]);
      else
        dy = 1.0;
      for (unsigned int j = 0; j < df.size (); j++)
        gsl_matrix_set (J, i, j, df[j]/dy);
    }
  }
  return GSL_SUCCESS;
}

#if GSL_MAJOR_VERSION > 1
#else
int astro::nlModeller::fJ
  (const gsl_vector * x, void * d, gsl_vector * f, gsl_matrix * J)
{
  mf (x, d, f);
  mJ (x, d, J);
  return GSL_SUCCESS;
}
#endif

astro::nlModeller::nlModeller (astro::data_model& data_in, bool jacobian):
data (data_in)
{
  fdf.f  = &mf;
  fdf.df = (jacobian) ? &mJ : NULL;
#if GSL_MAJOR_VERSION > 1
#else
  fdf.fdf = &fJ;
#endif
  fdf.n = data.n;
  fdf.p = data.p;
  fdf.params = &data;
}

#if GSL_MAJOR_VERSION > 1
void astro::nlModeller::operator ()
  (double * q, const unsigned int n_max, const double epsabs, const double epsrel)
{
  gsl_vector * q0 = gsl_vector_alloc (data.p);
  for (unsigned int i = 0; i < data.p; i++)
    gsl_vector_set (q0, i, q[i]);
  work = gsl_multifit_nlinear_alloc (type, &parameters, data.n, data.p);
  gsl_multifit_nlinear_init (q0, &fdf, work);
  int info = 0;
  int status = gsl_multifit_nlinear_driver
    (n_max, epsabs, epsrel, epsrel, NULL, NULL, &info, work);
  if (status == GSL_EMAXITER)
    throw std::runtime_error
      ("maximum number of iterations reached in nlModeller");
  for (unsigned int i = 0; i < data.p; i++)
    q[i] = gsl_vector_get (work->x, i);
}
#else
void astro::nlModeller::operator ()
  (double * q, const unsigned int n_max, const double epsabs, const double epsrel)
{
  gsl_vector_view p = gsl_vector_view_array (q, fdf.p);
  gsl_multifit_fdfsolver * solver =
    gsl_multifit_fdfsolver_alloc (type, fdf.n, fdf.p);
  gsl_multifit_fdfsolver_set (solver, &fdf, &p.vector);
  unsigned int iter = 0;
  int status = 0;
  do
  {
    iter++;
    gsl_multifit_fdfsolver_iterate (solver);
    status = gsl_multifit_test_delta (solver->dx, solver->x, epsabs, epsrel);
  }
  while (status == GSL_CONTINUE && iter < n_max);
  for (unsigned int i = 0; i < fdf.p; i++)
    q[i] = gsl_vector_get (solver->x, i);
}
#endif

void astro::nlModeller::get_cov_matrix
  (double * p, gsl_matrix * cov, const double epsrel)
{
  gsl_vector * x = gsl_vector_alloc (data.p);
  for (unsigned int i = 0; i < data.p; i++)
    gsl_vector_set (x, i, p[i]);
  gsl_matrix * J = gsl_matrix_alloc (data.n, data.p);
  mJ (x, &data, J);
#if GSL_MAJOR_VERSION > 1
  gsl_multifit_nlinear_covar (J, epsrel, cov);
#else
  gsl_multifit_covar (J, epsrel, cov);
#endif
  gsl_matrix_free (J);
  gsl_vector_free (x);
}
