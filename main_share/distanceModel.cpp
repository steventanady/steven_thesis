#include "distanceModel.h"

#include <astro/io/functionWriter.h>
#include "astro/utilities/utilities.h"
#include <astro/utilities/nlModeller.h>

distanceModel::distanceModel (dataSet * sample, quasarSample * qsample,
  astro::functionSet * basis_in, double exponent_in) : basis (basis_in)
{

  n_basis = basis->size ();
  a_min = basis->get_x_min ();
  a_max = basis->get_x_max ();
  exponent = exponent_in;

  n_data = 0.0;
  nsample = false;

  if (sample != NULL){

    n_data  = sample->size ();
    a_min_s = sample->get_a_min();
    nsample = true;

    scale    = sample->get_a ();
    distance = sample->get_distance ();
    name = sample->get_name ();

    Cov = sample->get_cov_matrix();

  }

  n_qso = 0;
  n_dataq = 0;
  quasars = false;

  if (qsample != NULL){

    n_dataq = qsample->size();
    if (nsample == false)
      n_qso = 2;
    else
      n_qso = 3;
    quasars = true;
    a_min_s = qsample->get_a_min();

    scaleq = qsample->get_a ();
    logf_x = qsample->get_logf_x ();
    logf_xerr = qsample->get_logf_xerr ();
    logf_uv = qsample->get_logf_uv ();
    logf_uverr = qsample->get_logf_uverr ();
    nameq = qsample->get_name ();

  }

  // do
  // {
    init_Coeffs ();
  // }
  // while (clean_Coeffs ());
}

distanceModel::~distanceModel ()
{
  gsl_matrix_free (cov);
}


void distanceModel::init_Fisher (data_model dm){

 J = gsl_matrix_alloc (dm.n,dm.p+dm.q);

 for (int i = 0; i < dm.n; i++){

   double result = 0.0;
   for (int j = 0; j < dm.p; j++){
     result += dm.basis->phi (dm.x[i], j)*p[j+dm.q];}
   double model=log10(result);

   if(dm.type[i]=="data/qso")
   {
     double delta_a = 1.0/dm.a_min-1.0;
     model += -log10(dm.a_min*dm.a_min*(1.0+dm.x[i]*delta_a)/sqrt(4.*M_PI)/3.0e9);
     // double r;
     // if (dm.nd == 0)
     //   r = (2.0-2.0*p[0])*model+dm.y[i]-p[0]*dm.z[i]+6.75;
     // else
     //   r = (2.0-2.0*p[0])*model+dm.y[i]-p[0]*dm.z[i]+1.0000*p[2];
     // double s_sq = pow(p[0]*dm.dz[i],2.0)+pow(dm.dy[i],2.0)+pow(p[1],2.0);

     gsl_matrix_set (J, i, 0, 2.*model+dm.z[i]);

     gsl_matrix_set (J, i, 1, 0.);

     if (nsample == true)
      gsl_matrix_set (J, i, 2, -1.0000);

     for (int m = 0; m < dm.p; m++){
       double t = -(2.-2.*p[0])*dm.basis->phi(dm.x[i], m)/result/log(10.);
       gsl_matrix_set(J, i, m+dm.q, t);}
   }
   else
   {
     if (quasars){
       gsl_matrix_set (J, i, 0, 0.);

       gsl_matrix_set (J, i, 1, 0.);

       if (nsample == true)
        gsl_matrix_set (J, i, 2, 0.);
     }

     for (int m = 0; m < dm.p; m++){
       double t = dm.basis->phi (dm.x[i], m)/result/log(10.);
       gsl_matrix_set(J, i, m+dm.q, t);}
   }
 }

 JTinvCov = gsl_matrix_alloc (dm.p+dm.q, dm.n);
 invtCov = gsl_matrix_calloc (dm.n, dm.n);

 for (int i = 0; i < dm.n; i++){
   if(dm.type[i]=="data/qso")
    gsl_matrix_set(invtCov,i,i,1./(pow(p[0]*dm.dz[i],2.0)+pow(dm.dy[i],2.0)+pow(p[1],2.0)));
   else{
     for (int j = 0; j < dm.nd; j++){
       gsl_matrix_set(invtCov,i,j,gsl_matrix_get(dm.invcov,i,j));
     }
   }
 }

 gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, J, invtCov, 0.0, JTinvCov);
 F = gsl_matrix_alloc (dm.p+dm.q, dm.p+dm.q);
 gsl_blas_dgemm
  (CblasNoTrans, CblasNoTrans, 1.0, JTinvCov, J, 0.0, F);

 gsl_matrix_free(J);
 gsl_matrix_free(JTinvCov);

  // Calculate Trace Part of the Fisher Matrix

  std::vector <gsl_matrix*> dCov;

  for (int i = 0; i < dm.p+dm.q; i++){
    gsl_matrix * a;
    a = gsl_matrix_alloc (dm.n, dm.n);
    gsl_matrix_set_all(a, 0.0);
    dCov.push_back(a);
  }

  for (int i = 0; i < dm.n; i++){

    if(dm.type[i]=="data/qso"){
     gsl_matrix_set(dCov[0],i,i,2.*p[0]*pow(dm.dz[i],2.0));
     gsl_matrix_set(dCov[1],i,i,2.*p[1]);}

  }

  std::vector <gsl_matrix*> invCovdCov;

  for (int i = 0; i < dm.p+dm.q; i++){
    gsl_matrix * a;
    a = gsl_matrix_alloc (dm.n, dm.n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invtCov, dCov[i], 0.0, a);
    invCovdCov.push_back(a);
  }

  gsl_matrix * Trace;
  Trace = gsl_matrix_alloc(dm.p+dm.q,dm.p+dm.q);

  for (int i = 0; i < dm.p+dm.q; i++){
    for (int j = 0; j < dm.p+dm.q; j++){
      gsl_matrix * temp;
      temp = gsl_matrix_alloc (dm.n, dm.n);
      gsl_blas_dgemm
        (CblasNoTrans, CblasNoTrans, 1.0, invCovdCov[i], invCovdCov[j], 0.0, temp);
      double trace = 0.0;
      for (int m = 0; m < dm.n; m++)
        trace += gsl_matrix_get(temp,m, m);
      gsl_matrix_set(Trace,i,j,trace);
      gsl_matrix_free(temp);
    }
  }

  gsl_matrix_scale(Trace,0.5);
  gsl_matrix_add(F,Trace);

  gsl_matrix_free(invtCov);
  gsl_matrix_free(Trace);
  for (int i = 0; i < dm.p+dm.q; i++){
    gsl_matrix_free(dCov[i]);
    gsl_matrix_free(invCovdCov[i]);
  }

}

double distanceModel::my_f (const gsl_vector *v, void *d){

  data_model dm = *reinterpret_cast<data_model *> (d);

  std::vector<double> p = astro::stdvec(v);

  double log_likelihood = 0.0;
  gsl_vector * mu = gsl_vector_alloc (dm.nd);

  for (int i = 0; i < dm.n; i++){

    double model = 0.0;
    for (int j = 0; j < dm.p; j++){
      model += dm.basis->phi (dm.x[i], j)*p[j+dm.q];}
    model = log10(model);

    if (dm.type[i]=="data/qso"){

      double delta_a = 1.0/dm.a_min-1.0;
      model += -log10(dm.a_min*dm.a_min*(1.0+dm.x[i]*delta_a)/sqrt(4.*M_PI)/3.0e9);
      double r;
      if (dm.nd == 0)
        r = (2.0-2.0*p[0])*model+dm.y[i]-p[0]*dm.z[i]+6.75;
      else
        r = (2.0-2.0*p[0])*model+dm.y[i]-p[0]*dm.z[i]+1.0000*p[2];
      double s_sq = pow(p[0]*dm.dz[i],2.0)+pow(dm.dy[i],2.0)+pow(p[1],2.0);
      log_likelihood += pow(r,2.)/s_sq+log(s_sq);

    }
    else
      gsl_vector_set(mu,i,model-dm.y[i]);

  }
  if(dm.nd!=0){
    gsl_vector * invcovTmu = gsl_vector_alloc (dm.nd);

    gsl_blas_dgemv (CblasTrans, 1.0, dm.invcov, mu, 0.0, invcovTmu);

    double temp;
    gsl_blas_ddot(invcovTmu, mu, &temp);

    gsl_vector_free (mu);
    gsl_vector_free (invcovTmu);

    log_likelihood += temp; // log(det(cov)) missing, but numerically unstable and just a constant
  }

  return log_likelihood;
}

void distanceModel::my_df (const gsl_vector *v, void *d, gsl_vector *df){

  data_model dm = *reinterpret_cast<data_model *> (d);

  std::vector<double> p = astro::stdvec(v);

  gsl_vector * mu = gsl_vector_alloc (dm.nd);
  std::vector <gsl_vector*> dmu;

  for (int i = 0; i < dm.p; i++){
    gsl_vector * a;
    a = gsl_vector_alloc (dm.nd);
    dmu.push_back(a);
  }

  for (int i = 0; i < dm.p+dm.q; i++)
    gsl_vector_set (df, i, 0.0);

  for (int i = 0; i < dm.n; i++){

    double result = 0.0;
    for (int j = 0; j < dm.p; j++){
      result += dm.basis->phi (dm.x[i], j)*p[j+dm.q];}
    double model=log10(result);

    if(dm.type[i]=="data/qso")
    {
      double delta_a = 1.0/dm.a_min-1.0;
      model += -log10(dm.a_min*dm.a_min*(1.0+dm.x[i]*delta_a)/sqrt(4.*M_PI)/3.0e9);
      double r;
      if (dm.nd == 0)
        r = (2.0-2.0*p[0])*model+dm.y[i]-p[0]*dm.z[i]+6.75;
      else
        r = (2.0-2.0*p[0])*model+dm.y[i]-p[0]*dm.z[i]+1.0000*p[2];
      double s_sq = pow(p[0]*dm.dz[i],2.0)+pow(dm.dy[i],2.0)+pow(p[1],2.0);

      gsl_vector_set(df,0,gsl_vector_get(df,0)-4.*r/s_sq*model+2.*p[0]*
      pow(dm.dz[i],2.)/s_sq*(1.-r*r/s_sq)-2.*r/s_sq*dm.z[i]);

      gsl_vector_set(df,1,gsl_vector_get(df,1)+2.*p[1]/s_sq*(1.-r*r/s_sq));

      if (dm.nd != 0)
        gsl_vector_set(df,2,gsl_vector_get(df,2)+2.*r/s_sq*1.0000);

      for (int m = 0; m < dm.p; m++){
        double t = 2.*r/s_sq*(2.-2.*p[0])*dm.basis->phi(dm.x[i], m)/result/log(10.);
        gsl_vector_set(df,m+dm.q,gsl_vector_get(df,m+dm.q)+t);}
    }
    else
    {
      gsl_vector_set(mu,i,model-dm.y[i]);

      for (int m = 0; m < dm.p; m++)
        gsl_vector_set(dmu[m],i,dm.basis->phi (dm.x[i], m)/result/log(10.));
    }
  }

    gsl_vector * invcovTmu = gsl_vector_alloc (dm.nd);
  if(dm.nd!=0){
    gsl_blas_dgemv (CblasTrans, 1.0, dm.invcov, mu, 0.0, invcovTmu);
  }

  for (int i = 0; i < dm.p; i++){
    double temp;
    gsl_blas_ddot(invcovTmu, dmu[i], &temp);
    gsl_vector_set(df,i+dm.q,gsl_vector_get(df,i+dm.q)+2.*temp);
  }
}

void::distanceModel::my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  my_df(x, params, df);
}

void distanceModel::init_Coeffs ()
{
  /*
  / Perform multidim minimization
  */

  if (n_data != 0){

  cov = gsl_matrix_alloc (n_data, n_data);
  gsl_matrix_memcpy(cov, Cov);

  for (int i = 0; i < n_data; i++){
    for (int j = 0; j < n_data; j++){
      gsl_matrix_set(cov,i,j,gsl_matrix_get(cov,i,j)
      /distance[i]/distance[j]/log(10.)/log(10.));
    }
  }
  invcov = gsl_matrix_alloc(n_data,n_data);

  gsl_matrix_memcpy (invcov, cov);
  gsl_linalg_cholesky_decomp1(invcov);
  gsl_linalg_cholesky_invert(invcov);

  // astro::inverse(cov, invcov);
  }

  bool cmb = false;
  bool pantheon = false;

  data_model d;
  d.n = n_data+n_dataq;
  d.nd = n_data;
  d.nq = n_qso;
  d.p = n_basis;
  d.q = n_qso;
  d.a_min = a_min_s;
  d.basis = basis;
  d.x = new double[n_data+n_dataq];
  d.y = new double[n_data+n_dataq];
  d.dy = new double [n_data+n_dataq];
  d.z = new double [n_data+n_dataq];
  d.dz = new double [n_data+n_dataq];
  d.cov = cov;
  d.invcov = invcov;
  d.type = new std::string [n_data+n_dataq];
  for (int i = 0; i < n_data; i++)
  {
    d.x[i] = scale[i];
    d.y[i] = log10(distance[i]);
    d.dy[i] = sqrt (gsl_matrix_get(Cov,i,i))/distance[i]/log(10.);
    d.z[i] = 0.0;
    d.dz[i] = 0.0;
    d.type[i] = name[i];
    if (name[i] == "data/cmbSample")
      cmb = true;
    if (name[i] == "data/pantheonSample")
      pantheon = true;
  }

  for (int i = 0; i < n_dataq; i++)
  {
    d.x[n_data+i] = scaleq[i];
    d.y[n_data+i] = logf_x[i];
    d.dy[n_data+i] = logf_xerr[i];
    d.z[n_data+i] = logf_uv[i];
    d.dz[n_data+i] = logf_uverr[i];
    d.type[n_data+i] = nameq[i];
  }

  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *m;

  //Initialize function for gsl multimin algorith
  gsl_multimin_function_fdf my_func;

  my_func.n = n_basis+n_qso;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = &d;

  // Starting values
  gsl_vector *x;
  x = gsl_vector_alloc (n_basis+n_qso);
  for (int i = 0; i < n_basis+n_qso; i++){
    gsl_vector_set(x,i,0.0001);
  }
  if (quasars) {
    gsl_vector_set (x, 0, 0.6);
    gsl_vector_set (x, 1, 0.2);
    if (nsample == true)
      gsl_vector_set (x, 2, 6.2);
    }

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  m = gsl_multimin_fdfminimizer_alloc (T, n_basis+n_qso);

  if(cmb)
    gsl_multimin_fdfminimizer_set (m, &my_func, x, 0.001, 1e-1);
  else if(pantheon)
    gsl_multimin_fdfminimizer_set (m, &my_func, x, 0.1, 1e-1);
  else
    gsl_multimin_fdfminimizer_set (m, &my_func, x, 0.01, 1e-1);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (m);

      if (status){
        std::cout<<status<<std::endl;
        break;
      }

      gsl_vector *test;
      test = gsl_vector_alloc(n_basis+n_qso);

      for (int i = 0; i < n_basis+n_qso; i++){
        gsl_vector_set(test, i, fabs(gsl_vector_get(m->gradient, i)
        *std::max(gsl_vector_get(m->x, i),1.))/std::max(fabs(m->f),1.));
      }

      if(cmb){
        if(gsl_vector_max(test)<1e-2)
          status = GSL_SUCCESS;
        else
          status = GSL_CONTINUE;}
      else if(pantheon){
        if(gsl_vector_max(test)<1e-4)
          status = GSL_SUCCESS;
        else
          status = GSL_CONTINUE;}
      else{
        if(gsl_vector_max(test)<1e-5)
          status = GSL_SUCCESS;
        else
          status = GSL_CONTINUE;}

      // status = gsl_multimin_test_gradient (test, 1e-4);
      std::cout << astro::str (m->x,5) << std::endl;

      if (status == GSL_SUCCESS){
        printf ("Minimum found at:\n");
        std::cout << astro::str (m->x,5) << std::endl;
        std::cout << "Likelihood: " << m->f << std::endl;
        std::cout << "Likelihood per dof: " << m->f/(n_data+n_dataq-n_basis+n_qso) << std::endl;
      }

    }
  while (status == GSL_CONTINUE && iter < 10000);


  p = astro::stdvec (m->x);

  init_Fisher (d);

  gsl_matrix * invF = gsl_matrix_alloc (n_basis+n_qso, n_basis+n_qso);

  astro::inverse (F, invF);

  /**
   * Copy results into standard vectors
   */
  dp = astro::stdvec (m->x);
  for (int i = 0; i < n_basis+n_qso; i++){
    dp[i] = sqrt(fabs(gsl_matrix_get(invF, i, i)));}
  std::cout << astro::str (dp) << std::endl;

  c.clear();
  for (int i = n_qso; i < n_basis+n_qso; i++){
    c.push_back(p[i]);}
  dc.clear();
  for (int i = n_qso; i < n_basis+n_qso; i++){
    dc.push_back(dp[i]);}

  gsl_multimin_fdfminimizer_free (m);
  gsl_vector_free (x);


  /**
   * Determine their uncertainties as follows: (1) diagonalize Fisher matrix;
   * (2) determine standard deviations in Fisher eigenframe; (3) rotate these
   * standard deviations back
   */

  v = gsl_vector_alloc (n_basis+n_qso);
  R = gsl_matrix_alloc (n_basis+n_qso, n_basis+n_qso);
  astro::eigen_system (F, v, R);


  gsl_vector * s = gsl_vector_alloc (n_basis+n_qso);
  for (int i = 0; i < n_basis+n_qso; i++)
    gsl_vector_set (s, i, 1./sqrt(gsl_vector_get (v, i)));
  gsl_vector * dp_temp = gsl_vector_alloc (n_basis+n_qso);
  gsl_blas_dgemv (CblasTrans, 1.0, R, s, 0.0, dp_temp);
  gsl_vector_free (s);

  dp = astro::stdvec (dp_temp);
  for (int i = 0; i < n_basis+n_qso; i++)
    dp[i] = fabs(dp[i]);
  dc.clear();
  for (int i = n_qso; i < n_basis+n_qso; i++){
    dc.push_back(dp[i]);}
  std::cout << astro::str (dp) << std::endl << std::endl;

  gsl_vector_free(v);
  gsl_vector_free(dp_temp);
  gsl_matrix_free(R);

  /* Marginalize over the nuisances parameters, needed by the quasars,
     by removing the corresponding rows and colummns from the inverse
     fisher matrix.
  */

  if (quasars){

    gsl_matrix * invFM;
    invFM = gsl_matrix_alloc (n_basis,n_basis);

    for (int i = 0; i < n_basis; i++){
      for (int j = 0; j < n_basis; j++){
        gsl_matrix_set(invFM,i,j,gsl_matrix_get(invF,i+n_qso,j+n_qso));
      }
    }

    gsl_vector * f;
    gsl_matrix * M;
    f = gsl_vector_alloc (n_basis);
    M = gsl_matrix_alloc (n_basis, n_basis);
    astro::eigen_system (invFM, f, M);

    gsl_vector * s = gsl_vector_alloc (n_basis);
    for (int i = 0; i < n_basis; i++)
      gsl_vector_set (s, i, sqrt(gsl_vector_get (f, i)));
    gsl_vector * dc_temp = gsl_vector_alloc (n_basis);
    gsl_blas_dgemv (CblasTrans, 1.0, M, s, 0.0, dc_temp);
    gsl_vector_free (s);

    dc.clear();
    for (int i = 0; i < n_basis; i++){
      dc.push_back(fabs(gsl_vector_get(dc_temp,i)));}
    std::cout << astro::str (dc) << std::endl;

    gsl_vector_free(f);
    gsl_vector_free(dc_temp);
    gsl_matrix_free(M);

  }
  //
  // delete p0;
  // gsl_matrix_free(invcov);



}

bool distanceModel::clean_Coeffs ()
{
  if (c.empty ())
    return false;
  bool flag = true;
  for (unsigned int i = 0; i < c.size (); i++)
    flag = flag && (fabs (c[i]) > dc[i]);
  if (!flag)
    n_basis--;
  return !flag;
}

double distanceModel::operator () (double a)
{
  double result = 0.0;
  for (int i = 0; i < n_basis; i++)
    result += c[i]*basis->phi (a, i);
  return log10(result);
}

double distanceModel::error (double a)
{
  double result = 0.0;
  for (int i = 0; i < n_basis; i++)
    result += dc[i]*basis->phi (a, i);
  return log10(result);
}

double distanceModel::prime (double a)
{
  double df = 0.0;

  for (int i = 0; i < n_basis; i++)
    df += basis->phi_prime (a, i)*c[i]/pow(operator()(a),10.0);

  return df;
}

double distanceModel::pprime (double a)
{
  double ddf = 0.0;

  for (int i = 0; i < n_basis; i++)
    ddf += basis->phi_pprime (a, i)*c[i]/pow(operator()(a),10.0)
           -pow(basis->phi_prime (a, i)*c[i]/pow(operator()(a),10.0),2.0);

  return ddf;
}

double distanceModel::get_c (int i)
{ return c[i]; }

std::vector<double> distanceModel::get_c ()
{ return c; }

double distanceModel::get_delta_c (int i)
{ return dc[i]; }

std::vector<double> distanceModel::get_delta_c ()
{ return dc; }

gsl_matrix * distanceModel::get_covariance_matrix ()
{ return cov; }

void distanceModel::print (std::string filename, int n)
{
  // std::ofstream output (filename);
  // astro::functionWriter write (&output);
  // write.push_back
  //   (std::bind (&distanceModel::operator (), this, std::placeholders::_1));
  // write.push_back
  //   (std::bind (&distanceModel::prime, this, std::placeholders::_1));
  // write (a_min, a_max);
  // output << std::endl << std::endl;
  // for (int i = 0; i < n_data; i++){
  //   output << scale[i] << "  " << distance[i] << "  ";
  //   output << sqrt(gsl_matrix_get(Cov,i,i)) << std::endl;}
  // output.close ();

  /*

  std::ofstream output2 ("distanceModel.d");
  astro::functionWriter write2 (&output2);
  write2.push_back
    (std::bind (&distanceModel::operator (), this, std::placeholders::_1));
  write2.push_back
    (std::bind (&distanceModel::prime, this, std::placeholders::_1));
  write2 (a_min, a_max);
  output2.close ();

  std::ofstream output3 ("distanceData.d");
  for (int i = 0; i < n_data; i++){
    output3 << scale[i] << "  " << log10(distance[i]) << "  ";
    output3 << sqrt(gsl_matrix_get(Cov,i,i))/distance[i]/log(10) << std::endl;}
  for (int i = 0; i < n_dataq; i++){
    output3 << scaleq[i] << "  " << 1./(2.-2.*p[0])*(p[0]*logf_uv[i]-logf_x[i])-p[0] << "  ";
    output3 << 0.0 << std::endl;}
  output3.close ();

  std::ofstream output4 ("residuals.d");
  for (int i = 0; i < n_data; i++){
    if (name[i]=="data/qso")
      output4 << scale[i] << "  " << log10(distance[i]-p[0])-operator ()(scale[i]) << std::endl;
    else
      output4 << scale[i] << "  " << log10(distance[i])-operator ()(scale[i]) << std::endl;}
  output4.close ();

  */

  double delta_a = 1.0/a_min_s-1.0;
  std::vector<double> f (n_data);
  double maxscale = a_max*(1.0-a_min_s)+a_min_s;
  double minscale = a_min*(1.0-a_min_s)+a_min_s;

  /*

  std::ofstream output5 ("luminosityDistance.d");
  for (int i = 0; i < n_data; i++){
    f[i]= a_min_s*a_min_s*(1.0+scale[i]*delta_a);
    output5 << scale[i]*(1.0-a_min_s)+a_min_s << "  " << (distance[i])/f[i] << "  ";
    output5 << sqrt(gsl_matrix_get(Cov,i,i))/f[i] << std::endl;}
  for (int i = 0; i < n_dataq; i++) {
    if (nsample == true)
      output5 << scaleq[i]*(1.0-a_min_s)+a_min_s << "  " << pow(10.,1./(2.-2.*p[0])*(p[0]*logf_uv[i]-logf_x[i]-p[2]*1.0000))/sqrt(4.*M_PI)/3.0e9 << "  ";
    else
      output5 << scaleq[i]*(1.0-a_min_s)+a_min_s << "  " << pow(10.,1./(2.-2.*p[0])*(p[0]*logf_uv[i]-logf_x[i]-6.75))/sqrt(4.*M_PI)/3.0e9 << "  ";
    output5 << 0.0 << std::endl;}
  output5.close ();

  */

  std::ofstream output6 ("luminosityDistanceModel.d");
  astro::functionWriter write6 (&output6);
  write6.push_back ([&] (double a) { return pow(10.0,operator ()((a-a_min_s)/(1.0-a_min_s)))*pow(a_min_s,-2.)/(1.0+((a-a_min_s)/(1.0-a_min_s))*delta_a)*3.0e9; });
  write6.push_back ([&] (double a) { return pow(10.0,error ((a-a_min_s)/(1.0-a_min_s)))*pow(a_min_s,-2.)/(1.0+((a-a_min_s)/(1.0-a_min_s))*delta_a)*3.0e9; });
  write6 (minscale, maxscale, 1028, astro::LOG_SPACING);
  output6.close ();

  if (quasars){
    std::ofstream output7 ("qsoSample.d");
    for (int i = 0; i < n_dataq; i++) {
      output7 << 1./(scaleq[i]*(1.0-a_min_s)+a_min_s)-1. << "  ";
      if (nsample == true)
        output7 << 5.*(1./(2.-2.*p[0])*(p[0]*logf_uv[i]-logf_x[i]-p[2]*1.0000)-0.5*log10(4.*M_PI))-5. << "  ";
      else
        output7 << 5.*(1./(2.-2.*p[0])*(p[0]*logf_uv[i]-logf_x[i]-6.75)-0.5*log10(4.*M_PI))-5. << "  ";
      output7 << 5./(2.-2.*p[0])*sqrt(pow(((1./(1.-p[0])*(p[0]*logf_uv[i]-logf_x[i])+
                    logf_uv[i])*dp[0]),2.)+pow((p[0]*logf_uverr[i]),2.)+pow(logf_xerr[i],2.))<< " ";
      output7 << 5./(2.-2.*p[0])*sqrt(pow(((1./(1.-p[0])*(p[0]*logf_uv[i]-logf_x[i])+
                    logf_uv[i])*dp[0]),2.)+pow((p[0]*logf_uverr[i]),2.)+pow(logf_xerr[i],2.)+pow(p[1],2.))<< std::endl;}
    output7.close ();
  }


}

void distanceModel::report
  (std::string filename, std::string input, double a_min, double exponent)
{
  std::ofstream output (filename);
  output
    << "Report card for data sample: " << input << std::endl
    << std::endl
    << "a_min: " << a_min << std::endl
    << "prefactor exponent: " << exponent << std::endl
    << std::endl
    << "coefficients:  " << astro::str (c) << std::endl
    << "uncertainties: " << astro::str (dc) << std::endl
    << std::endl;
    if (quasars){
      output << "gamma: " << p[0] << " +- " << dp[0] << std::endl
      << "delta: " << p[1] << " +- " << dp[1] << std::endl;
      if (nsample == true)
        output << "beta:  " << p[2] << " +- " << dp[2] << std::endl;
      output << std::endl;
    }
  output.close ();
}
