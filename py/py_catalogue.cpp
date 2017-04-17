//
// wrapping lib/catalogue.cpp
//
#include "msg.h"
#include "catalogue.h"
#include "error.h"
#include "py_catalogue.h"
#include "py_util.h"
#include "py_assert.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#ifdef DOUBLEPRECISION
#define NPY_FLOAT_TYPE NPY_DOUBLE
#else
#define NPY_FLOAT_TYPE NPY_FLOAT
#endif

PyMODINIT_FUNC
py_catalogue_module_init()
{
  import_array();

  return NULL;
}


static void py_catalogue_free(PyObject *obj);
static void py_catalogue_file_free(PyObject* obj);

//
// wrapping Catalogue class
//
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

//
// Wrapping CatalogueFile class
//
PyObject* py_catalogue_file_ascii_alloc(PyObject* self, PyObject* args)
{
  // Create a new CatalogueFile object
  // Args:
  //     filename
  //     ipos: a sequence of column indeces for positions
  //     iweights: a sequence of column indeces for weights
  //     inbar: index for nbar
  //     Pest:  Pest for FKP weight

  PyObject *bytes;
  PyObject *py_ipos, *py_iweights;
  int inbar;
  double Pest;
		     
  if(!PyArg_ParseTuple(args, "O&OOid", PyUnicode_FSConverter, &bytes,
			 &py_ipos, &py_iweights, &inbar, &Pest)) {
    return NULL;
  }

  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  std::vector<int> ipos, iweights;
  sequence_to_vector(py_ipos, ipos);
  sequence_to_vector(py_iweights, iweights);

  CatalogueFile* const f=
    new CatalogueFileAscii(filename, ipos, iweights, inbar, Pest);
  
  Py_DECREF(bytes);

  return PyCapsule_New(f, "_CatalogueFile", py_catalogue_file_free);
}

void py_catalogue_file_free(PyObject* obj)
{
  CatalogueFile* const f=
    (CatalogueFile*) PyCapsule_GetPointer(obj, "_CatalogueFile");
  py_assert_void(f);

  delete f;
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

PyObject* py_catalogue_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_cat;

  if(!PyArg_ParseTuple(args, "O", &py_cat)) {
    return NULL;
  }

  Catalogue const * const cat=
    (Catalogue const *) PyCapsule_GetPointer(py_cat, "_Catalogue");
  py_assert_ptr(cat);

  const int nd=2;
  const int ncol= sizeof(Particle)/sizeof(Float);
  npy_intp dims[]= {(npy_intp) cat->size(), ncol};

  return PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT_TYPE,
				   (Float*) cat->data());
}


