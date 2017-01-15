//
// wrapping lib/catalogue.cpp
//
#include "msg.h"
#include "catalogue.h"
#include "error.h"
#include "py_catalogue.h"
#include "py_assert.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

PyMODINIT_FUNC
py_catalogue_module_init()
{
  msg_printf(msg_debug, "Module py_catalogue initialised for numpy array.\n");
  import_array();

  return NULL;
}


static void py_catalogue_free(PyObject *obj);

PyObject* py_catalogue_alloc(PyObject* self, PyObject* args)
{
  // _catalogue_alloc()
  // Create a new catalogue object
  Catalogue* const cat= new Catalogue();

  return PyCapsule_New(cat, "_Catalogue", py_catalogue_free);
}

void py_catalogue_free(PyObject *obj)
{
  // Delete the catalogue object
  // Called automatically by Python
  Catalogue* const cat= (Catalogue*) PyCapsule_GetPointer(obj, "_Catalogue");
  py_assert_void(cat);

  delete cat;
}

PyObject* py_catalogue_read_text(PyObject* self, PyObject* args)
{
  // _catalogue_read_ascii(_cat,filename)
  // Read particles from an ascii file into the catalogue
  PyObject *py_cat, *bytes;

  if(!PyArg_ParseTuple(args, "OO&", &py_cat, PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }

  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  Catalogue* const cat= (Catalogue*) PyCapsule_GetPointer(py_cat, "_Catalogue");
  py_assert_ptr(cat);

  try {
    catalogue_read_text(cat, filename);
  }
  catch(IOError) {
    Py_DECREF(bytes);
    PyErr_SetString(PyExc_IOError, "Catalogue file error");
    return NULL;
  }

  Py_DECREF(bytes);

  Py_RETURN_NONE;
}

PyObject* py_catalogue_len(PyObject* self, PyObject* args)
{
  // _catalogue_len(_cat)
  // Return the number of particles in the catalogue
  PyObject *py_cat;

  if(!PyArg_ParseTuple(args, "O", &py_cat)) {
    return NULL;
  }

  Catalogue const * const cat=
    (Catalogue const *) PyCapsule_GetPointer(py_cat, "_Catalogue");
  py_assert_ptr(cat);


  return Py_BuildValue("n", cat->size());
}


