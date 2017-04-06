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
    if(std::is_same<T, int>::value)
      v.push_back(PyLong_AsLong(elem));
    else if(std::is_same<T, double>::value)
      v.push_back(PyFloat_AsDouble(elem));
    else {
      PyErr_SetString(PyExc_TypeError,
		      "The value in the sequence is not the type expected");
    }
  }
}

	
    

#endif
