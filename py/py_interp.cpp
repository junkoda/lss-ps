#include <cassert>
#include "py_interp.h"

using namespace std;


Interp::Interp(const vector<double>& v_x, const vector<double>& v_y)
{
  assert(v_x.size() > 0); // for front() and back() to exist
  
  const size_t n= v_x.size();
  assert(v_y.size() == n);
  
  acc= gsl_interp_accel_alloc();
  spline= gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, v_x.data(), v_y.data(), n);

  x_min = v_x.front();
  x_max = v_x.back();
}

Interp::~Interp()
{
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

double Interp::operator()(const double x)
{
  if(x < x_min || x > x_max)
    return 0.0;
  
  return gsl_spline_eval(spline, x, acc);
}

double Interp::f(const double x)
{
  if(x < x_min || x > x_max)
    return 0.0;

  return gsl_spline_eval(spline, x, acc);
}
