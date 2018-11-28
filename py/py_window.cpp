#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "py_util.h"
#include "py_window.h"
#include "py_interp.h"
#include "py_assert.h"
#include "config.h"

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

/*
template<typename float_type>
double compute_statistics_from_array(float_type const * const xyz,
				     const size_t xyz_stride,
				     float_type const * const weight,
				     const size_t weight_stride,
				     const size_t np,
				     const int n)
*/


PyObject* py_window_compute_shotnoise(PyObject* self, PyObject* args)
{
  PyObject *py_xyz, *py_weight, *py_nbar;
  double x0[3];
  int n;
  if(!PyArg_ParseTuple(args, "iOOOddd",
		       &n,
		       &py_xyz, &py_weight, &py_nbar,
		       x0, x0 + 1, x0 + 2))
    return NULL;

  // xyz array
  Py_buffer buf_xyz;
  if(PyObject_GetBuffer(py_xyz, &buf_xyz, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;
  
  if(buf_xyz.ndim != 2) {
    PyErr_SetString(PyExc_TypeError, "Expected a 2-dimensional array for xyz");
    return NULL;
  }

  // weight array
  Py_buffer buf_weight;
  Py_ssize_t zero[]= {(Py_ssize_t) 0};
  
  if(py_weight == Py_None) {
    buf_weight.buf= nullptr;
    buf_weight.suboffsets= nullptr;
    buf_weight.strides= zero;
  }
  else {
    if(PyObject_GetBuffer(py_weight, &buf_weight, PyBUF_FORMAT | PyBUF_FULL_RO)
       == -1) return NULL;

    if(buf_weight.ndim != 1) {
      PyErr_SetString(PyExc_TypeError,
		      "Expected a 1-dimensional array for weight");
      return NULL;
    }
  }

  // nbar
  Py_buffer buf_nbar;
  if(py_nbar == Py_None) {
    buf_nbar.buf= nullptr;
    buf_nbar.suboffsets= nullptr;
    buf_nbar.strides= zero;
  }
  else {
    if(PyObject_GetBuffer(py_nbar, &buf_nbar,
			  PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
      return NULL;

    if(buf_nbar.ndim != 1) {
      PyErr_SetString(PyExc_TypeError,
		      "Expected a 1-dimensional array for nbar");
      return NULL;
    }
  }

  if(buf_xyz.suboffsets || buf_weight.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "_window_compute_shotnoise cannot handle array with suboffsets");
    return NULL;
  }

  if(buf_weight.buf && (buf_weight.shape[0] != buf_xyz.shape[0])) {
    PyErr_SetString(PyExc_TypeError, "Length of xyz and weight arrays differ");
    return NULL;
  }

  if(strcmp(buf_xyz.format, "d") == 0) {
    if(buf_xyz.strides[1] != sizeof(double)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }
  }

  double* weight= (double*) buf_weight.buf;
  double* xyz= (double*) buf_xyz.buf;
  double* nbar= (double*) buf_nbar.buf;

  const size_t xyz_stride= buf_xyz.strides[0];
  const size_t weight_stride= buf_weight.strides[0];
  const size_t nbar_stride= buf_nbar.strides[0];
  
  size_t np= static_cast<size_t>(buf_xyz.shape[0]);

  long double w2_sum= 0.0, nw2_sum= 0.0;

  for(size_t i=0; i<np; ++i) {
    Float x[3];
    double w = weight == nullptr ? 1.0 : *weight;
      
    x[0] = xyz[0] - x0[0];
    x[1] = xyz[1] - x0[1];
    x[2] = xyz[2] - x0[2];
    
    Float r= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    const double w2 = w*w;
    w2_sum += w2/pow(r, n);
    nw2_sum += (*nbar)*w2;
    
    xyz      = (double*) ((char*) xyz    + xyz_stride);
    if(weight)
      weight = (double*) ((char*) weight + weight_stride);
    if(nbar)
      nbar   = (double*) ((char*) nbar   + nbar_stride);
  }

  double shotnoise= static_cast<double>(w2_sum/nw2_sum);
  return Py_BuildValue("d", shotnoise);
}

