#include <cstdio>
#include <cassert>
#include "py_cola_binary.h"
#include "py_assert.h"

namespace {
template<typename float_type>
void load_cola_particles(FILE* fp, float_type* p, size_t np, size_t stride)
{
  float header[6];
  int ret= fread(header, sizeof(float), 6, fp); assert(ret == 6);

  int np_header;
  ret= fread(&np_header, sizeof(int), 1, fp); assert(ret == 1);
  assert(np == static_cast<size_t>(np_header));
  
  float x[6];
  for(size_t i=0; i<np; ++i) {
    ret= fread(x, sizeof(float), 6, fp); assert(ret == 6);
    for(int k=0; k<6; ++k)
      p[k]= x[k]; //static_cast<float_type>(x[k]);
    
    p = (float_type*) ((char*) p + stride);
  }

  int np_check;
  ret= fread(&np_check, sizeof(int), 1, fp); assert(ret == 1);
  assert(np == static_cast<size_t>(np_check));
}
}

PyObject* py_cola_binary_load_particles_header(PyObject* self, PyObject* args)
{
  PyObject* bytes;
  
  if(!PyArg_ParseTuple(args, "O&",
		       PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }

  // Decode filename
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  FILE* fp= fopen(filename, "r");
  Py_DECREF(bytes);

  if(fp == 0) {
    PyErr_SetNone(PyExc_FileNotFoundError);
    return NULL;
  }

  int np;
  float header[6];
  
  fread(header, sizeof(float), 6, fp);
  fread(&np, sizeof(int), 1, fp);

  fclose(fp);

  return Py_BuildValue("(idddddd)", np,
		       header[0], header[1], header[2],
		       header[3], header[4], header[5]);
}
  
PyObject* py_cola_binary_load_particles(PyObject* self, PyObject* args)
{
  // _cola_binary_load_particles(filename, a)

  PyObject *py_a;
  PyObject* bytes;

  if(!PyArg_ParseTuple(args, "O&O",
		       PyUnicode_FSConverter, &bytes, &py_a)) {
    return NULL;
  }

  // Decode filename
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  FILE* fp= fopen(filename, "r");
  Py_DECREF(bytes);

  if(fp == 0) {
    PyErr_SetNone(PyExc_FileNotFoundError);
    return NULL;
  }
  
  // Decode array
  Py_buffer a;
  if(PyObject_GetBuffer(py_a, &a, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;
  
  if(a.ndim != 2) {
    PyErr_SetString(PyExc_TypeError, "Expected a 2-dimensional array for xyz");
    return NULL;
  }

  if(a.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "cannot handle array with suboffsets");
    return NULL;
  }

  py_assert_ptr(a.shape[1] == 6);

  if(strcmp(a.format, "d") == 0) {
    py_assert_ptr(a.strides[1] == sizeof(double));
    load_cola_particles<double>(fp, (double*) a.buf, a.shape[0], a.strides[0]);
  }
  else if(strcmp(a.format, "f") == 0) {
    py_assert_ptr(a.strides[1] == sizeof(float));
    load_cola_particles<float>(fp, (float*) a.buf, a.shape[0], a.strides[0]);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "Expected an array of floats or doubles");
    return NULL;
  }

  PyBuffer_Release(&a);


  Py_RETURN_NONE;
}
