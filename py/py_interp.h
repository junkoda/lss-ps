#ifndef PY_INTERP_PY
#define PY_INTREP_PY 1

#include <vector>
#include <gsl/gsl_spline.h>
#include "Python.h"

class Interp {
 public:
  Interp(const std::vector<double>& v_x, const std::vector<double>& v_y);
  ~Interp();
  double operator()(const double x);
  double f(const double x);
  double x_min, x_max;
 private:
  gsl_spline* spline;
  gsl_interp_accel* acc;
};

#endif
