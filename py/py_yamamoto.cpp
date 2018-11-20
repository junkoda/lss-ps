#include <vector>

#include "config.h" // Float
#include "moment.h" // Moment1, Moment2, ...
#include "yamamoto.h"
#include "py_yamamoto.h"
#include "py_util.h"

using namespace std;

PyObject* py_yamamoto_compute_moment_x(PyObject* self, PyObject* args)
{
  // _yamamoto_compute_moment_x(_grid_delta, ij, _grid_moment)
  // grid_delta
  // ij (list): list of indeces len(ij) = 2 or 4
  //            e.g. Q_ij(x) for [i, j], Q_ijkl(x) for [i, j, k, l]

  PyObject *py_grid_delta, *py_list, *py_x0, *py_grid_moment;
  int n;
  
  if(!PyArg_ParseTuple(args, "OOOOi",
		       &py_grid_delta, &py_x0, &py_list, &py_grid_moment,
		       &n))
    return NULL;

  Grid const * const delta=
    (Grid const *) PyCapsule_GetPointer(py_grid_delta, "_Grid");
  py_assert_ptr(delta);

  Grid* const moment=
    (Grid*) PyCapsule_GetPointer(py_grid_moment, "_Grid");
  py_assert_ptr(moment);

  vector<int> v;
  sequence_to_vector<int>(py_list, v);

  vector<Float> x0;
  sequence_to_vector<Float>(py_x0, x0);

  if(v.size() == 1) {
    if(n == 0)
      moment_compute_x(delta, x0.data(), Moment1(v[0]), moment);
    else
      moment_compute_x(delta, x0.data(), WindowMoment1(v[0], n), moment);
  }
  else if(v.size() == 2) {
    if(n == 0)
      moment_compute_x(delta, x0.data(), Moment2(v[0], v[1]), moment);
    else
      moment_compute_x(delta, x0.data(), WindowMoment2(v[0], v[1], n), moment);
  }
  else if(v.size() == 3) {
    if(n == 0)
      moment_compute_x(delta, x0.data(), Moment3(v[0], v[1], v[2]), moment);
    else
      moment_compute_x(delta, x0.data(), WindowMoment3(v[0], v[1], v[2], n),
		       moment);
      
  }
  else if(v.size() == 4) {
    if(n == 0)
      moment_compute_x(delta, x0.data(), Moment4(v[0], v[1], v[2], v[3]),
		       moment);
    else
      moment_compute_x(delta, x0.data(),
		       WindowMoment4(v[0], v[1], v[2], v[3], n), moment);
  }
  else {
    PyErr_SetString(PyExc_ValueError, "length of indecies is nither 2 or 4");
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

  if(v.size() == 1) {
    moment_compute_k(delta, Moment1(v[0]), coef, moment);
  }
  else if(v.size() == 2) {
    moment_compute_k(delta, Moment2(v[0], v[1]), coef, moment);
  }
  else if(v.size() == 3) {
    moment_compute_k(delta, Moment3(v[0], v[1], v[2]), coef, moment);
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




