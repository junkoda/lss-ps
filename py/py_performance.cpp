//
// Perform basic operations to measure basic performance
//

#include <iostream>
#include "msg.h"
#include "py_assert.h"
#include "py_util.h"

using namespace std;

template<typename float_type>
void sum_from_array(float_type const * xyz,
		    const size_t xyz_stride,
		    const size_t np,
		    double sums[])
{
  assert(np > 0);

  double sum_x=0.0, sum_y=0.0, sum_z=0.0;

#ifdef _OPENMP
  #pragma omp parallel for reduction(+ : sum_x, sum_y, sum_z)
#endif
  for(size_t i=0; i<np; ++i) {
    float_type const * const x=
      (float_type const * const)((char const * const) xyz + xyz_stride*i);

    sum_x += x[0];
    sum_y += x[1];
    sum_z += x[2];
  }

  sums[0]= sum_x;
  sums[1]= sum_y;
  sums[2]= sum_z;
}


PyObject* py_performance_xyz_sum(PyObject* self, PyObject* args)
{
  // _min_from_array(xyz)
  // Compute minimum of xyz
  // Args:
  //     xyz:    The array of xyz
  //
  PyObject *py_xyz;

  if(!PyArg_ParseTuple(args, "O", &py_xyz))
    return NULL;

  //
  // Decode array information
  //

  // xyz array
  Py_buffer xyz;
  if(PyObject_GetBuffer(py_xyz, &xyz, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;

  if(xyz.ndim != 2) {
    PyErr_SetString(PyExc_TypeError, "Expected a 2-dimensional array for xyz");
    return NULL;
  }

  if(xyz.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "_mass_assignment_from_array cannot handle array with suboffsets");
    return NULL;
  }

  if(xyz.strides[0] == 0)
    Py_RETURN_NONE;

  double sums[3];

  if(strcmp(xyz.format, "d") == 0) {
    if(xyz.strides[1] != sizeof(double)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    sum_from_array<double>((double*) xyz.buf, xyz.strides[0], xyz.shape[0],
                           sums);
  }
  else if(strcmp(xyz.format, "f") == 0) {
    if(xyz.strides[1] != sizeof(float)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    sum_from_array<float>((float*) xyz.buf, xyz.strides[0], xyz.shape[0],
			  sums);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "Expected an array of floats or doubles");
    return NULL;
  }


  PyBuffer_Release(&xyz);

  return Py_BuildValue("(ddd)", sums[0], sums[1], sums[2]);
}

