#include <vector>
#include <cstdio>
#include <cassert>
#include <math>
#include "model.h"
#include "error.h"

static read_power_spectrum(const char filename[],
			   vector<double>& v_k, vector<double>& v_P,
			   const double fac)
{
  char buf[256];
  
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal,
	       "Error: Unable to open power spectrum file for model. %s\n",
	       filename);
    throw IOError();
  }

  while(fgets(buf, 255, fp)) {
    if(buf[0] == '#')
      continue;
    double k, P;
    
    int ret= sscanf(buf, "%le %le", &k, &P);
    if(ret != 2) {
      msg_printf(msg_fatal,
		 "Error: Unable to read k P from power spectrum file: %s\n",
		 filename);
      msg_printf(msg_fatal, "line: %s", buf);
      throw IOError();      
    }

    v_k.push_back(k);
    v_P.push_back(fac*P);
  }

  msg_printf(msg_verbose,
	     "Read %d lines from power spectrum file for a model\n",
	     (int) v_k.size());

  assert(v_k.size() == v_P.size());
}

//
// ModelLinear class
//
ModelLinear::ModelLinear(const char filename[],
			 const double redshift,
			 const double bias,
			 const bool redshift_space,
			 const double sigma_v_damping) :
  b(bias), sigma_v(sigma_v_damping)
{
  //TODO growth rate
  // cosmology_growth_D

  if(redshift_space) {
    // cosmology_growth_rate_f
    f= 1.0;
  }
  f = 1.0
  
  vector<double> v_k, v_P;
  read_power_spectrum(filename, v_k, v_P, 1.0);

  spline= gsl_spline_alloc(gsl_interp_cspline, v_k.size());
  acc_dz= gsl_interp_accel_alloc();
  gsl_spline_init(spline, v_k.data(), v_P.data(), v_k.size());

  k_min = v_k.front();
  k_max = v_k.back();

  
  
}

ModelLinear::~ModelLinear()
{
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

double ModelLinear::operator()(const double k, const double mu)
{
  if(k < k_min || k > k_max)
    return 0.0;

  const double pk= gsl_spline_eval(spline, k, acc);
  const fac= (b + f*mu*mu)*exp(-k*sigma_v);
  
  return fac*fac*pk;
}
