#include "py_util.h"
#include "error.h"

using namespace std;

static void decode_array(const char name[],
			 PyObject* py_obj, Py_buffer* buf,
			 const Py_ssize_t len=0, const Py_ssize_t ncol=1,
			 const bool read_only=true)
{
  // Common procedure for using Python array
  //
  // name: name of the array for error message
  // py_obs: array object
  // buf: resulting buffer object
  // len: expected length; raise error if the length of the array is not len
  
  char msg[128];

  int flag = 0;
  if(read_only)
    flag= PyBUF_FULL_RO;
  else
    flag= PyBUF_FULL;
  
  if(PyObject_GetBuffer(py_obj, buf, PyBUF_FORMAT | flag) == -1)
      throw TypeError();

  if(ncol == 1) {
    if(buf->ndim != 1) {
      sprintf(msg, "Expected a 1-dimensional array for %s", name);
      PyErr_SetString(PyExc_TypeError, msg);
      throw TypeError();
    }
  }
  else if(ncol > 1) {
    if(buf->ndim != 2) {
      sprintf(msg, "Expected a 2-dimensional array for %s", name);
      PyErr_SetString(PyExc_TypeError, msg);
      throw TypeError();
    }
  }


  if(strcmp(buf->format, "d") != 0) {
    sprintf(msg, "Expected an array of double for %s: %s", name, buf->format);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  if(len > 0 && buf->shape[0] != len) {
    sprintf(msg, "Expected the length arrays of %d for %s: %d",
	    (int) len, name, (int) buf->shape[0]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }
  
  if(ncol > 0 && buf->shape[1] != ncol) {
    sprintf(msg, "Expected number of columns %d for %s: %d",
	    (int) ncol, name, (int) buf->shape[2]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

}


void py_util_array_as_vector(const char name[],
			     PyObject* py_array,
			     vector<double>& v,
			     const Py_ssize_t len_expect)
{
  // 1-dimensional array of double as a vector
  Py_buffer buf;

  decode_array(name, py_array, &buf, len_expect, 0, true);

  double const * x= (double const *) buf.buf;
  const size_t n= buf.shape[0];
  const size_t stride= buf.strides[0];

  v.reserve(n);
  
  for(size_t i=0; i<n; ++i) {
    v.push_back(*x);
    x = (double const *) ((char const *) x + stride);
  }

  PyBuffer_Release(&buf);
}
	
    
