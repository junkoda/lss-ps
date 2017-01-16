//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_config.h"
#include "py_catalogue.h"
#include "py_grid.h"
#include "py_mass_assignment.h"

using namespace std;

static PyMethodDef methods[] = {
  {"_config_sizeof_float", py_config_sizeof_float, METH_VARARGS,
   "return sizeof(Float)"},

  {"_catalogue_alloc", py_catalogue_alloc, METH_VARARGS,
   "allocate a new _cat opbject"},
  {"_catalogue_len", py_catalogue_len, METH_VARARGS,
   "return number of particles in the catalogue"},
  {"_catalogue_read_text", py_catalogue_read_text, METH_VARARGS,
   "_catalogue_read_text(_cat, filename) read particle data"},
  {"_catalogue_asarray", py_catalogue_asarray, METH_VARARGS,
   "_catalogue_asarrary(_cat) return catalogue as a np.array"},

  {"_grid_alloc", py_grid_alloc, METH_VARARGS,
   "_grid_alloc(nc) allocate a new grid object"},
  {"_grid_fft", py_grid_fft, METH_VARARGS,
   "_grid_fft(_grid); apply FFT to the grid"},
  {"_grid_nc", py_grid_nc, METH_VARARGS, "return nc"},
  {"_grid_mode", py_grid_mode, METH_VARARGS, "return grid mode"},
  {"_grid_fx_asarray", py_grid_fx_asarray, METH_VARARGS,
   "return real-space grid as an np.array"},
  {"_grid_fk_asarray", py_grid_fk_asarray, METH_VARARGS,
   "return Fourier-space grid as an np.array"},

  {"_mass_assignment_cic", py_mass_assignment_cic, METH_VARARGS,
   "assign mass to a grid using CIC"},

  //{py_power_spectrum_alloc(PyObject* self, PyObject* args)
  PyObject* py_power_spectrum_alloc(PyObject* self, PyObject* args)
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_lssps", // name of this module
  "A package for power spectrum cosmology:large-scale structure", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__lssps(void) {
  py_catalogue_module_init();
  py_grid_module_init();
    
  
  return PyModule_Create(&module);
}
