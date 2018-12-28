#include <iostream> // DEBUG
#include "py_buffer.h"

//
// Wrap Buffer protocol
//

using namespace std;

Buffer::Buffer(char const * const name_, PyObject* const py_array) :
  name(name_)
{
  // may throw TypeError
  
  char msg[64];

  if(py_array == Py_None) {
    pybuf.buf= nullptr;
    pybuf.ndim= 0;
    str_format= "";
  }
  else {
    if(PyObject_GetBuffer(py_array, &pybuf, PyBUF_FORMAT | PyBUF_FULL) == -1) {
      sprintf(msg, "Expected a buffer protocol for %.16s", name);
      PyErr_SetString(PyExc_TypeError, msg);
      throw TypeError();
    }

    const char c= pybuf.format[0];

    // first character may indicate the byteorder
    if(c == '@' || c == '=' || c == '<') { // assuming native is little endian
      str_format= pybuf.format[1];
    }
    else {
      str_format= pybuf.format;
    }
  }


}

Buffer::~Buffer()
{
  if(!pybuf.buf)
    PyBuffer_Release(&pybuf);
}

Py_ssize_t Buffer::shape(const int i) const
{
  // may throw ValueError
  //cerr << "shape() " << str_format << endl;//DEBUG
  //cerr << "ndim " << pybuf.ndim << endl;
  //cerr << "shape " << pybuf.shape[0] << endl;

  char msg[96];
  if(i < 0 || i >= pybuf.ndim) {
    sprintf(msg, "shape(%d) undefined for a %d-dimensional array for %.16s",
	    i, pybuf.ndim, name);
    PyErr_SetString(PyExc_IndexError, msg);
    throw ValueError();
  }

  return pybuf.shape[i];
}

void Buffer::assert_shape(const int n_dim,
		     const Py_ssize_t len, const Py_ssize_t ncol)
{
  // assert that the buffer has the expected shape
  char msg[64];

  if(n_dim == 1) {
    if(pybuf.ndim != 1) {
      sprintf(msg, "Expected a 1-dimensional array for %.16s", name);
      PyErr_SetString(PyExc_AssertionError, msg);
      throw AssertionError();
    }

    if(pybuf.shape[0] != len) {
      sprintf(msg, "Expected a length %d for %.16s", (int) len, name);
      PyErr_SetString(PyExc_AssertionError, msg);
      throw AssertionError();
    }
  }
  else {
    if(pybuf.ndim != 2) {
      sprintf(msg, "Expected a 2-dimensional array for %.16s", name);
      PyErr_SetString(PyExc_AssertionError, msg);
      throw AssertionError();
    }
    
    if(len > 0 && pybuf.shape[0] != len) {
      sprintf(msg, "Expected %d rows for %.16s: %d",
	      (int) len, name, (int) pybuf.shape[0]);
      PyErr_SetString(PyExc_AssertionError, msg);
      throw AssertionError();
    }
    else if(ncol > 0 && pybuf.shape[1] != ncol) {
      sprintf(msg, "Expected %d columns for %.16s: %d",
	      (int) ncol, name, (int) pybuf.shape[1]);
      PyErr_SetString(PyExc_AssertionError, msg);
      throw AssertionError();
    }
  }
}

void Buffer::assert_format(char const * const fmt) {
  char msg[64];
  if(strcmp(pybuf.format, fmt) != 0) {
    sprintf(msg, "Expected an array of float for %.16s: %.4s",
	    name, pybuf.format);
    PyErr_SetString(PyExc_AssertionError, msg);
    throw AssertionError();
  }
}
  
