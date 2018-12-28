#include "py_buffer.h"
#include "py_util.h"
#include "error.h"

using namespace std;

void py_util_array_as_vector(const char name[],
			     PyObject* py_array,
			     vector<double>& v,
			     const Py_ssize_t len_expect)
{
  // 1-dimensional array of double as a vector
  // Exceptions:
  //   may throw AssersionError
  Buffer buf(name, py_array);

  if(len_expect > 0)
    buf.assert_shape(1, len_expect);

  double const * x= buf.buf<double>(); //(double const *) buf.buf;
  const size_t n= buf.shape(0);
  const size_t stride= buf.stride(0);

  v.reserve(n);
  
  for(size_t i=0; i<n; ++i) {
    v.push_back(*x);
    x = (double const *) ((char const *) x + stride);
  }
}
	
    
void py_util_vector_as_array(const char name[], const vector<double>& v,
			     PyObject* py_array)
{
  // copy vector v content to array py_array
  // The length of the array must be the same as that of the vector
  //
  // Exceptions:
  //   may throw TypeError()

  Buffer buf(name, py_array);
  buf.assert_shape(1, v.size());

  double* x= buf.buf<double>();
  const size_t n= buf.shape(0);
  const size_t stride= buf.stride(0);

  for(size_t i=0; i<n; ++i) {
    *x= v[i];
    x= (double *) ((char const *) x + stride);
  }
}
