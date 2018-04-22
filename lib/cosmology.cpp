#include "cosmology.h"

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

static const double cH0inv= 2997.92458; // c/H0 [h^-1 Mpc]

static double omega_m= 0.0;
static double z_max, d_max= 0.0;

static gsl_interp_accel *acc_dz= 0;
static gsl_interp_accel *acc_zd= 0;
static gsl_spline* spline_dz= 0;
static gsl_spline* spline_zd= 0;

static double growth_int(double a, void* param);
static double growth(const double a, const double omega_m);


void cosmology_init(const double omega_m_, const double z_max_, const size_t n)
{
  omega_m= omega_m_;
  z_max= z_max_;

  assert(omega_m > 0.0);

  // Precompute a table of comoving distances
  double* const d= (double*) malloc(sizeof(double)*n*2); assert(d);
  double* const z= d + n;

  d[0]= 0.0; z[0]= 0.0;

  for(size_t i=1; i<n; ++i) {
    z[i]= z_max*i/(n-1);
    double a= 1.0/(1.0 + z[i]);
    
    d[i]= cosmology_compute_comoving_distance(a);
  }
  d_max= d[n-1];

  spline_dz= gsl_spline_alloc(gsl_interp_cspline, n);
  acc_dz= gsl_interp_accel_alloc();
  gsl_spline_init(spline_dz, d, z, n);

  spline_zd= gsl_spline_alloc(gsl_interp_cspline, n);
  acc_zd= gsl_interp_accel_alloc();
  gsl_spline_init(spline_zd, z, d, n);

  free(d);
}

void cosmology_free()
{
  if(omega_m > 0.0) {
    gsl_spline_free(spline_dz);
    gsl_spline_free(spline_zd);
    gsl_interp_accel_free(acc_dz);
    gsl_interp_accel_free(acc_zd);
  }
}

  
double cosmology_omega_m()
{
  return omega_m;
}

static double distance_integrand(double a, void* params) {
  // 1/[ a^2 H(a)/H0 ]
  const double om= *(double*) params;
  return 1.0/sqrt(om*a + (1.0 - om)*(a*a*a*a));
}

double cosmology_compute_comoving_distance(const double a)
{
  // Light emitted at comoving distance r at scale facgor a reach
  // the observer now
  //
  // r = \int c*dt/aH(a) = \int c*da/[a^2 H(a)]
  //   = c/H0 \int da/[a^2 sqrt(omega_m*a + omega_lambda*a^4)]
  
  const int nwork= 1000;
  gsl_integration_workspace* w= gsl_integration_workspace_alloc(nwork);

  gsl_function F;
  F.function= &distance_integrand;
  F.params= (void*) &omega_m;

  double result, error;
  gsl_integration_qags(&F, a, 1.0, 0, 1.0e-8, nwork, w, &result, &error);

  gsl_integration_workspace_free(w);

  return cH0inv*result;
}

double cosmology_distance_redshift(const double d)
{
  assert(omega_m > 0.0);
  if(d > d_max)
    return -1;
  
  return gsl_spline_eval(spline_dz, d, acc_dz);
}

double cosmology_redshift_distance(const double z)
{
  assert(omega_m > 0.0);
  assert(z <= z_max);
  
  return gsl_spline_eval(spline_zd, z, acc_zd);
}

//
// growth factor
//
double cosmology_growth_D(const double a, const double omega_m)
{
  // Growth factor D(a)
  return growth(a, omega_m)/growth(1.0, omega_m);
}

double cosmology_growth_rate_f(const double a, const double omega_m)
{
  // Growth rate f=dlnD/dlna
  const double D0= growth(1.0, omega_m);
  const double hubble_a= sqrt(omega_m/(a*a*a) + 1.0 - omega_m);
  const double d= growth(a, omega_m)/D0;

  const double f_ex= 1.0/(d*D0*a*a*hubble_a*hubble_a)
                     - 1.5*omega_m/(hubble_a*hubble_a*a*a*a);

  return f_ex;
}


double growth_int(double a, void* param)
{
  // integrand for growth factor
  const double om= *(double*)param;
  return pow(a/(om + (1.0 - om)*a*a*a), 1.5);
}


double growth(const double a, const double omega_m)
{
  // Compute integral of D_tilda(a) -- unnormalised growth factor
  // growth_factor = D_tilda(a)/D_tilda(1.0)

  const double hubble_a= sqrt(omega_m/(a*a*a) + 1.0 - omega_m);

  const int worksize= 1000;
  double result, abserr;

  gsl_integration_workspace *workspace=
    gsl_integration_workspace_alloc(worksize);

  gsl_function F;
  F.function = &growth_int;
  F.params = (void*) &omega_m;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, worksize, GSL_INTEG_GAUSS41,
		      workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}
