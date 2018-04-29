#include <cstdio>
#include "msg.h"
#include "py_util.h"
#include "py_assert.h"
#include "py_gadget_file.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

struct gadget_file_header{
  int      np[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int np_total[6];
  int      flag_cooling;
  int      num_files;
  double   boxsize;
  double   omega0;
  double   omega_lambda;
  double   hubble_param; 
  int flag_stellarage;
  int flag_metals;
  unsigned int np_total_highword[6];
  int  flag_entropy_instead_u;
  char fill[60];
  //char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
           /* fills to 256 Bytes */
};

PyMODINIT_FUNC
py_gadget_file_module_init()
{
  import_array();

  return NULL;
}


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

  gadget_file_header h;
  ret= fread(&h, sizeof(gadget_file_header), 1, fp);
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

  gadget_file_header h;
  ret= fread(&h, sizeof(gadget_file_header), 1, fp);
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
