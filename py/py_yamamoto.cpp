#include <vector>

#include "moment.h"
#include "yamamoto.h"
#include "py_yamamoto.h"
#include "py_util.h"

using namespace std;

PyObject* py_yamamoto_compute_moment_x(PyObject* self, PyObject* args)
{
  // _yamamoto_compute_moment_x(_grid_delta, ij, _grid_moment)
  // grid_delta
  // ij (list): list of indeces len(ij) = 2 or 4

  PyObject *py_grid_delta, *py_list, *py_grid_moment;

  if(!PyArg_ParseTuple(args, "OOO",
		       &py_grid_delta, &py_list, &py_grid_moment))
    return NULL;

  Grid const * const delta=
    (Grid const *) PyCapsule_GetPointer(py_grid_delta, "_Grid");
  py_assert_ptr(delta);

  Grid* const moment=
    (Grid*) PyCapsule_GetPointer(py_grid_moment, "_Grid");
  py_assert_ptr(moment);

  vector<int> v;
  sequence_to_vector<int>(py_list, v);

  if(v.size() == 2) {
    moment_compute_x(delta, Moment2(v[0], v[1]), moment);
  }
  else if(v.size() == 4) {
    moment_compute_x(delta, Moment4(v[0], v[1], v[2], v[3]), moment);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "length of indecies is nither 2 or 4");
    return NULL;
  }

  Py_RETURN_NONE;
}

PyObject* py_yamamoto_compute_moment_k(PyObject* self, PyObject* args)
{
  // _yamamoto_compute_moment_x(_grid_delta, ij, _grid_moment)
  // grid_delta
  // ij (list): list of indeces len(ij) = 2 or 4

  PyObject *py_grid_delta, *py_list, *py_grid_moment;
  double coef;
  
  if(!PyArg_ParseTuple(args, "OOdO",
		       &py_grid_delta, &py_list, &coef,
		       &py_grid_moment))
    return NULL;

  Grid const * const delta=
    (Grid const *) PyCapsule_GetPointer(py_grid_delta, "_Grid");
  py_assert_ptr(delta);

  Grid* const moment=
    (Grid*) PyCapsule_GetPointer(py_grid_moment, "_Grid");
  py_assert_ptr(moment);

  vector<int> v;
  sequence_to_vector<int>(py_list, v);

  if(v.size() == 2) {
    moment_compute_k(delta, Moment2(v[0], v[1]), coef, moment);
  }
  else if(v.size() == 4) {
    moment_compute_k(delta, Moment4(v[0], v[1], v[2], v[3]), coef, moment);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "length of indecies is nither 2 or 4");
    return NULL;
  }

  Py_RETURN_NONE;
}




