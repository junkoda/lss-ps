#ifndef MODEL_H
#define MODEL_H 1

#include <gsl/gsl_spline.h>

  
class ModelLinear
{
 public:
  ModelLinear(const char filename[],
	      const double redshift,
	      const double bias,
	      const double f_growth_rate,
	      const double sigma_v_damping);
  ~ModelLinear();
  double operator()(const double k, const double mu);

  double b, f, sigma_v;
  double k_min, k_max;
 private:
  gsl_interp_accel* acc;
  gsl_spline* spline;
};


 
#endif
