#include <cstdio>
#include "msg.h"
#include "error.h"
#include "py_util.h"
#include "py_assert.h"
#include "gadget_file.h"
#include "py_gadget_file.h"

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include "numpy/arrayobject.h"

static void py_gadget_file_free(PyObject *obj);

PyObject* py_gadget_file_read_np(PyObject* self, PyObject* args)
{
  //_gadget_file_read_np(filename)
  PyObject *bytes;
  if(!PyArg_ParseTuple(args, "O&",
		       PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }
  
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  FILE* const fp= fopen(filename, "r");
  
  if(fp == 0) {
    msg_printf(msg_error, "Error: Gadget file not found: %s", filename);
    PyErr_SetString(PyExc_FileNotFoundError, "File not found");
    Py_DECREF(bytes);
    return NULL;
  }

  int sep_begin, sep_end;
  int ret= fread(&sep_begin, sizeof(int), 1, fp); py_assert_ptr(ret == 1);
  if(sep_begin != 256) {
    msg_printf(msg_error, "Error: Unexpeced format for Gadget binary file: %s",
	       filename);
    PyErr_SetString(PyExc_ValueError, "Expected value in gadget binary file");
    Py_DECREF(bytes);
    return NULL;
  }

  GadgetFileHeader h;
  ret= fread(&h, sizeof(GadgetFileHeader), 1, fp);
  py_assert_ptr(ret == 1);
  
  ret= fread(&sep_end, sizeof(int), 1, fp);
  py_assert_ptr(ret == 1);
  py_assert_ptr(sep_begin == sep_end);

  fclose(fp);
  Py_DECREF(bytes);

  int* np= h.np;
  return Py_BuildValue("iiiiii",
		       np[0], np[1], np[2], np[3], np[4], np[5]);
}

PyObject* py_gadget_file_read_header(PyObject* self, PyObject* args)
{
  //_gadget_file_read_np(filename)
  PyObject *bytes;
  if(!PyArg_ParseTuple(args, "O&",
		       PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }
  
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  FILE* const fp= fopen(filename, "r");
  
  if(fp == 0) {
    msg_printf(msg_error, "Error: Gadget file not found: %s", filename);
    PyErr_SetString(PyExc_FileNotFoundError, "File not found");
    Py_DECREF(bytes);
    return NULL;
  }

  int sep_begin, sep_end;
  int ret= fread(&sep_begin, sizeof(int), 1, fp); py_assert_ptr(ret == 1);
  if(sep_begin != 256) {
    msg_printf(msg_error, "Error: Unexpeced format for Gadget binary file: %s",
	       filename);
    PyErr_SetString(PyExc_ValueError, "Expected value in gadget binary file");
    Py_DECREF(bytes);
    return NULL;
  }

  GadgetFileHeader h;
  ret= fread(&h, sizeof(GadgetFileHeader), 1, fp);
  py_assert_ptr(ret == 1);
  
  ret= fread(&sep_end, sizeof(int), 1, fp);
  py_assert_ptr(ret == 1);
  py_assert_ptr(sep_begin == sep_end);

  fclose(fp);
  Py_DECREF(bytes);

  unsigned long long npt[6];
  for(int i=0; i<6; ++i) {
    npt[i]=
      (static_cast<unsigned long long>(h.np_total_highword[i]) << 32) +
       static_cast<unsigned long long>(h.np_total[i]);
  }

  int* np= h.np;
  double* m= h.mass;
  return Py_BuildValue("(iiiiii)(KKKKKK)(dddddd)ddidddd",
		       np[0], np[1], np[2], np[3], np[4], np[5],
		       npt[0], npt[1], np[2], npt[3], npt[4], npt[5],
		       m[0], m[1], m[2], m[3], m[4], m[5],
		       h.time, h.redshift, h.num_files,
		       h.boxsize, h.omega0, h.omega_lambda, h.hubble_param);
}


PyObject* py_gadget_file_alloc(PyObject* self, PyObject* args)
{
  //_gadget_file_alloc(filename)
  PyObject *bytes;
  if(!PyArg_ParseTuple(args, "O&",
		       PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }
  
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  GadgetFile* const gf= new GadgetFile(filename);
  
  return PyCapsule_New(gf, "_GadgetFile", py_gadget_file_free);
}

void py_gadget_file_free(PyObject *obj)
{
  // Delete the power spectrum object
  // Called automatically by Python
  GadgetFile* const gf=
    (GadgetFile*) PyCapsule_GetPointer(obj, "_GadgetFile");
  py_assert_void(gf);

  delete gf;
}

PyObject* py_gadget_file_open(PyObject* self, PyObject* args)
{
  PyObject *py_gf;

  if(!PyArg_ParseTuple(args, "O", &py_gf)) {
    return NULL;
  }

  GadgetFile* const gf=
    (GadgetFile*) PyCapsule_GetPointer(py_gf, "_GadgetFile");
  py_assert_ptr(gf);

  try {
    gf->open();
  }
  catch(FileNotFoundError) {
    PyErr_SetString(PyExc_FileNotFoundError, "File not found");
    return NULL;
  }    

  Py_RETURN_NONE;
}

PyObject* py_gadget_file_close(PyObject* self, PyObject* args)
{
  PyObject *py_gf;

  if(!PyArg_ParseTuple(args, "O", &py_gf)) {
    return NULL;
  }

  GadgetFile* const gf=
    (GadgetFile*) PyCapsule_GetPointer(py_gf, "_GadgetFile");
  py_assert_ptr(gf);

  gf->close();

  Py_RETURN_NONE;
}

PyObject* py_gadget_file_read(PyObject* self, PyObject* args)
{
  //
  // _gadget_file_read(_gf, component, ibegin, iend, a)
  // Args:
  //     _gf (_GadgetFile): _GadgetFile pointer
  //     component (str): one-character string 'x' or 'v'
  //     ibegin - iend: index range of reading particles
  //     a (array): output array of components
  PyObject *py_gf, *py_array;
  char* component;
  int ibegin, iend;

  if(!PyArg_ParseTuple(args, "OsiiO", &py_gf, &component,
		       &ibegin, &iend,  &py_array)) {
    return NULL;
  }

  GadgetFile* const gf=
    (GadgetFile*) PyCapsule_GetPointer(py_gf, "_GadgetFile");
  py_assert_ptr(gf);

  Py_buffer a;
  if(PyObject_GetBuffer(py_array, &a, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;

  if(a.ndim != 2) {
    PyErr_SetString(PyExc_TypeError, "Expected a 2-dimensional array for _gadget_file_read");
    return NULL;
  }

  if(strcmp(a.format, "f") != 0) {
    PyErr_SetString(PyExc_TypeError, "Expected a float array for _gadget_file_read");
    return NULL;
  }

  gf->read(component[0], ibegin, iend, a.strides[0], (float*) a.buf);

  Py_RETURN_NONE;
}
