#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "py_util.h"
#include "py_window.h"
#include "py_interp.h"
#include "py_assert.h"

using namespace std;

namespace {
struct Params {
  int l;
  double r;
  Interp* interp;
};
}

double integrand_bessel(double k, void* vparams)
{
  // k^2 j_l(kr) f(k)
  Params* params= (Params*) vparams;
  double kr= k*params->r;
  return k*k*gsl_sf_bessel_jl(params->l, kr)*params->interp->f(k);
}

PyObject* py_window_cquag_bessel_transform(PyObject* self, PyObject* args)
{
  int l;
  PyObject *py_r, *py_k, *py_f, *py_result;
  if(!PyArg_ParseTuple(args, "OOOiO", &py_r, &py_k, &py_f, &l, &py_result))
    return NULL;

  vector<double> v_r, v_k, v_f, v_result;

  py_util_array_as_vector("r", py_r, v_r);
  py_util_array_as_vector("k", py_k, v_k); 
  py_util_array_as_vector("f", py_f, v_f, v_k.size());
  py_assert_ptr(v_k.size() > 0);

  v_result.resize(v_r.size(), 0.0);


  Interp interp(v_k, v_f);

  gsl_integration_cquad_workspace* const work=
    gsl_integration_cquad_workspace_alloc(100);

  Params params;
  params.interp= &interp;
  params.l= l;

  gsl_function F;
  F.function = &integrand_bessel;
  F.params = &params;

  const double k_min= v_k.front();
  const double k_max= v_k.back();
  const double fac= 1.0/(2.0*M_PI*M_PI);
  
  for(size_t i=0; i<v_r.size(); ++i) {
    params.r= v_r[i];
    
    double result, err;
    gsl_integration_cquad(&F, k_min, k_max, 0.0, 1.0e-5, work,
			  &result,&err, NULL);

    v_result[i]= fac*result;
  }

  gsl_integration_cquad_workspace_free(work);

  py_util_vector_as_array("result", v_result, py_result);
  
  Py_RETURN_NONE;
}



