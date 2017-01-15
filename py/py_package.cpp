//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_config.h"
#include "py_catalogue.h"

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
  
  return PyModule_Create(&module);
}
