#ifndef _PY_BUFFER_H
#define _PY_BUFFER_H 1

#include <string>
#include <vector>
#include <typeinfo>
#include "Python.h"
#include "error.h"

	       
class Buffer {
 public:
  Buffer(char const * const name_, PyObject* const py_array);
  ~Buffer();

  Buffer(const Buffer&) = delete;
  Buffer& operator=(const Buffer&) = delete;

  int ndim() const {
    return pybuf.ndim;
  }
  
  char const * format() const {
    return pybuf.format;
  }

  Py_ssize_t const * strides() const {
    return pybuf.strides;
  }

  Py_ssize_t shape(const int i) const;

  

  template <typename T> T* buf() const {
    char msg[64];
    
    if(str_format == "d") {
      if(typeid(T) != typeid(double)) {
	sprintf(msg, "Expected an array of double for %.16s but it is %.4s",
		name, pybuf.format);
	PyErr_SetString(PyExc_TypeError, msg);
	throw TypeError();
      }
    }
    else if(str_format == "f") {
      if(typeid(T) != typeid(float)) {
	sprintf(msg, "Expected an array of float for %.16s but it is %.4s",
		name, pybuf.format);
	PyErr_SetString(PyExc_TypeError, msg);
	throw TypeError();
      }
    }  
    else {
      sprintf(msg, "Array %.16s has an unknown dtype %.4s",
	      name, str_format.c_str());
      PyErr_SetString(PyExc_TypeError, msg);
      throw TypeError();
    }
    
    return reinterpret_cast<T*>(pybuf.buf);
  }

  /*
    template<bool Cond>
    auto f(){
    if constexpr(Cond)
    return 0;
    else
    return 3.14f;
    }
  */
  /*
  template<> constexpr float* buf() const {
    return reinterprete_cast<float*>(pybuf.buf);
  }
  template<> constexpr double* buf() const {
    return reinterprete_cast<double*>(pybuf.buf);
  }
  */
  
  void assert_shape(const int n_dim,
		    const Py_ssize_t len, const Py_ssize_t ncol=0);

  void assert_format(char const * const fmt);
 private:
  Py_buffer pybuf;
  char const * const name;
  std::string str_format;
  // format string
  // d for double, f for float, ... 
  // https://docs.python.org/3/library/struct.html#module-struct
 };

#endif
