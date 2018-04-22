#include "model.h"
#include "py_assert.h"
#include "py_model.h"
#include "py_power_spectrum.h"

static void py_model_free(PyObject *obj);

void py_model_free(PyObject *obj)
{
  // Delete the model object
  // Called automatically by Python
  Model* const ps_model=
    (Model*) PyCapsule_GetPointer(obj, "_Model");
  py_assert_void(ps_model);

  delete ps_model;
}

PyObject* py_model_linear_alloc(PyObject* self, PyObject* args)
{
  // _model_linear_alloc(filename, omega_m, z, b, redshift_space, sigma_v)
  PyObject* bytes;
  double omega_m, z, b, sigma_v;
  int redshift_space;

  if(!PyArg_ParseTuple(args, "O&dddid",
		       PyUnicode_FSConverter, &bytes,
		       &omega_m, &z, &b, &redshift_space, &sigma_v)) {
    return NULL;
  }
  
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);
  
  Model* ps= new ModelLinear(filename, omega_m, z, b, redshift_space, sigma_v);

  Py_DECREF(bytes);

  return PyCapsule_New(ps, "_Model", py_model_free);
}

PyObject* py_model_call(PyObject* self, PyObject* args)
{
  // _model_call(_model, k, mu)
  // return P_model(k, mu)
  
  PyObject* py_model;
  double k, mu;
  if(!PyArg_ParseTuple(args, "Odd", &py_model, &k, &mu)) {
    return NULL;
  }

  Model* const ps=
    (Model*) PyCapsule_GetPointer(py_model, "_Model");
  py_assert_ptr(ps);

  return Py_BuildValue("d", (*ps)(k, mu));
}

PyObject* py_model_compute_discrete_multipoles(PyObject* self, PyObject* args)
{
  // _model_compute_discrite_multipoles(_model, nc, boxsize, k_min, k_max, dk,
  // mas_window_deg
  // Evaluate the model on 3D grid, and then compute multipoles
  PyObject* py_model;
  int nc;
  double boxsize, k_min, k_max, dk;
  int mas_window_nc, mas_window_deg;

  if(!PyArg_ParseTuple(args, "Oiddddii",
		       &py_model, &nc, &boxsize, &k_min, &k_max, &dk,
		       &mas_window_nc, &mas_window_deg)) {
    return NULL;
  }
  
  Model const * const model=
    (Model const*) PyCapsule_GetPointer(py_model, "_Model");
  py_assert_ptr(model);

  PowerSpectrum* const ps=
    model_compute_discrete_multipoles(nc, boxsize,
				      k_min, k_max, dk,
				      *model, mas_window_nc, mas_window_deg);

  
  return PyCapsule_New(ps, "_PowerSpectrum", py_power_spectrum_free);
}
