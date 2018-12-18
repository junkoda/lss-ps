#ifndef PY_UTIL_H
#define PY_UTIL_H 1

#include <vector>
#include <type_traits>
#include "Python.h"
#include "py_assert.h"

template<typename T>
void sequence_to_vector(PyObject* seq, std::vector<T>& v)
{
  //
  // Add the contens of the Python squence to C++ vector v
  // T is either int or double
  if(PySequence_Check(seq) == false) {
    PyErr_SetString(PyExc_TypeError,
		    "sequence is not given to sequene_to_vector");
  }

  size_t n = (size_t) PySequence_Length(seq);
  v.reserve(v.size() + n);

  for(size_t i=0; i<n; ++i) {
    PyObject* elem = PySequence_GetItem(seq, i); py_assert_void(elem);
    // new reference
    if(std::is_same<T, int>::value)
      v.push_back(PyLong_AsLong(elem));
    else if(std::is_same<T, double>::value)
      v.push_back(PyFloat_AsDouble(elem));
    else {
      PyErr_SetString(PyExc_TypeError,
		      "The value in the sequence is not the type expected");
    }
    Py_DECREF(elem);
  }
}


void py_util_array_as_vector(const char name[],
			     PyObject* py_array,
			     std::vector<double>& v,
			     const Py_ssize_t len_expect=0);

void py_util_vector_as_array(const char name[], const std::vector<double>& v,
			     PyObject* py_array);

void py_util_decode_array(const char name[],
			  PyObject* py_obj, Py_buffer* buf,
			  const Py_ssize_t len=0, const Py_ssize_t ncol=1,
			  const bool read_only=true);

#endif
