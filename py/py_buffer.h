#ifndef _PY_BUFFER_H
#define _PY_BUFFER_H 1

#include <vector>
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

  template <typename T>
  T* buf() const {
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
};

#endif
