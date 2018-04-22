#ifndef MODEL_H
#define MODEL_H 1

#include <vector>
#include <gsl/gsl_spline.h>
#include "power_spectrum.h"
#include "grid.h"

//
// Model classes
//
class Model
{
 public:
  virtual ~Model() = default;
  virtual double operator()(const double k, const double mu) const = 0;
};

class ModelLinear : public Model
{
 public:
  ModelLinear(const char filename[],
	      const double omega_m,
	      const double redshift,
	      const double bias,
	      const bool redshift_space,
	      const double sigma_v_damping);
  ModelLinear(const ModelLinear&) = delete;
  ModelLinear& operator=(const ModelLinear &) = delete;

  virtual ~ModelLinear();
  virtual double operator()(const double k, const double mu) const;

  double b, f, sigma_v;
  double k_min, k_max;
 private:
  gsl_interp_accel* acc;
  gsl_spline* spline;
};

PowerSpectrum* model_compute_discrete_multipoles(const int nc,
						 const Float boxsize,
						 const double k_min,
						 const double k_max,
						 const double dk,
						 const Model& f,
						 const int mas_window_nc,
						 const int mas_window_deg);

double model_apply_window_3d(const Model& f,
			     Grid const * const grid,
			     const double k);

 
#endif
