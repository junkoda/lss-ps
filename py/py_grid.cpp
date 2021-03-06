//
// wrapping lib/grid.cpp
//
#include <iostream>
#include <cmath>

#include "msg.h"
#include "config.h"
#include "grid.h"
#include "error.h"
#include "py_assert.h"
#include "py_util.h"
#include "py_interp.h"
#include "py_grid.h"

using namespace std;

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#ifdef DOUBLEPRECISION
#define NPY_FLOAT_TYPE NPY_DOUBLE
#define NPY_COMPLEX_TYPE NPY_COMPLEX128
#else
#define NPY_FLOAT_TYPE NPY_FLOAT
#define NPY_COMPLEX_TYPE NPY_COMPLEX64
#endif

PyMODINIT_FUNC
py_grid_module_init()
{
  import_array();

  return NULL;
}


static void py_grid_free(PyObject *obj);

template<typename T>
static void write_float_to_binary(FILE* fp, Float const * p, const size_t nc)
{
  const size_t ncz= 2*(nc/2+1);
  vector<T> v;
  v.reserve(nc);

  for(size_t i=0; i<nc*nc; ++i) {
    v.clear();
    for(size_t j=0; j<nc; ++j)
      v.push_back(static_cast<T>(p[j]));
    int ret= fwrite(v.data(), sizeof(T), nc, fp);
    assert(ret == static_cast<int>(nc));
    p += ncz;
  }  
}

//
// class Grid
//

PyObject* py_grid_alloc(PyObject* self, PyObject* args)
{
  // _grid_alloc(nc)
  // Create a new grid object

  int nc;
  if(!PyArg_ParseTuple(args, "i", &nc)) {
    return NULL;
  }

  //msg_printf(msg_verbose, "Allocating a grid nc = %d\n", nc);

  Grid* grid= 0;
  try {
    grid= new Grid(nc);
  }
  catch(MemoryError) {
    PyErr_SetString(PyExc_MemoryError, "Grid memory error");
    return NULL;
  }

  return PyCapsule_New(grid, "_Grid", py_grid_free);
}

void py_grid_free(PyObject *obj)
{
  // Delete the grid object; called automatically by Python
  Grid* const grid= (Grid*) PyCapsule_GetPointer(obj, "_Grid");
  py_assert_void(grid);

  delete grid;
}

PyObject* py_grid_fft(PyObject* self, PyObject* args)
{
  // Fourier transform the grid
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  try {
    grid->fft();
  }
  catch(AssertionError) {
    PyErr_SetString(PyExc_MemoryError, "Grid FFT error");
  }

  Py_RETURN_NONE;
}

PyObject* py_grid_fft_inverse(PyObject* self, PyObject* args)
{
  // Inverse Fourier transform the grid
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  try {
    grid->fft_inverse();
  }
  catch(AssertionError) {
    PyErr_SetString(PyExc_MemoryError, "Grid Inverse FFT error");
  }

  Py_RETURN_NONE;
}

PyObject* py_grid_nc(PyObject* self, PyObject* args)
{
  // _grid_nc(_grid)
  // Return the number of grid points per dimension
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("n", grid->nc);
}

PyObject* py_grid_get_mode(PyObject* self, PyObject* args)
{
  // _grid_mode(_grid)
  // Return the mode of the grid real-space / fourier-space
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  if(grid->mode == GridMode::real_space)
    return Py_BuildValue("s", "real-space");
  else if(grid->mode == GridMode::fourier_space)
    return Py_BuildValue("s", "fourier-space");

  Py_RETURN_NONE;
}

PyObject* py_grid_set_mode(PyObject* self, PyObject* args)
{
  // _grid_mode(_grid, mode)
  // set mode
  // mode: 0 unknown
  //       1 real-space
  //       2 fourier-space
  PyObject *py_grid;
  int mode;

  if(!PyArg_ParseTuple(args, "Oi", &py_grid, &mode)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->mode= static_cast<GridMode>(mode);

  Py_RETURN_NONE;
}


PyObject* py_grid_get_boxsize(PyObject* self, PyObject* args)
{
  // _grid_get_boxsize(_grid)
  // Return boxsize
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("d", grid->boxsize);
}

PyObject* py_grid_set_boxsize(PyObject* self, PyObject* args)
{
  // _grid_set_boxsize(_grid, boxsize)
  // Set boxsize
  PyObject *py_grid;
  double boxsize;

  if(!PyArg_ParseTuple(args, "Od", &py_grid, &boxsize)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->boxsize= boxsize;

  Py_RETURN_NONE;
}

PyObject* py_grid_get_x0(PyObject* self, PyObject* args)
{
  // _grid_get_x0(_grid, x0)
  // Return x0_box
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("(d,d,d)",
		       grid->x0_box[0], grid->x0_box[1], grid->x0_box[2]);  
}

PyObject* py_grid_set_x0(PyObject* self, PyObject* args)
{
  // _grid_set_x0(_grid, x0, y0, z0)
  // Set x0_box
  PyObject *py_grid;
  double x0[3];

  if(!PyArg_ParseTuple(args, "Oddd", &py_grid, x0, x0+1, x0+2)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  for(int k=0; k<3; ++k)
    grid->x0_box[k]= x0[k];

  Py_RETURN_NONE;
}

PyObject* py_grid_get_offset(PyObject* self, PyObject* args)
{
  // _grid_get_offset(_grid)
  // Return offset
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("d", grid->offset);
}

PyObject* py_grid_set_offset(PyObject* self, PyObject* args)
{
  // _grid_set_offset(_grid, offset)
  // set offset
  PyObject *py_grid;
  double offset;

  if(!PyArg_ParseTuple(args, "Od", &py_grid, &offset)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->offset= offset;

  Py_RETURN_NONE;
}

PyObject* py_grid_get_sums(PyObject* self, PyObject* args)
{
  // _grid_get_sums(_grid)
  // Return sum of weights, total_weight, w2_sum, nw2_sum
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("(dddk)",
		       (double) grid->total_weight,
		       (double) grid->w2_sum,
		       (double) grid->nw2_sum,
		       (unsigned long) grid->np);
}

PyObject* py_grid_set_sums(PyObject* self, PyObject* args)
{
  // _grid_set_sum(_grid, total_weight, w2_sum, nw2_sum, np)
  PyObject *py_grid;
  double w_sum, w2_sum, nw2_sum;
  unsigned long np;

  if(!PyArg_ParseTuple(args, "Odddk", &py_grid,
		       &w_sum, &w2_sum, &nw2_sum, &np)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->total_weight = w_sum;
  grid->w2_sum = w2_sum;
  grid->nw2_sum = nw2_sum;
  grid->np = np;

  Py_RETURN_NONE;
}

PyObject* py_grid_get_nmas(PyObject* self, PyObject* args)
{
  // _grid_get_nmas(_grid)
  // Return the degree of mass assighment
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("i", grid->n_mas);
}

PyObject* py_grid_set_nmas(PyObject* self, PyObject* args)
{
  // _grid_set_nmas(_grid, nmas)
  PyObject *py_grid;
  int n_mas;

  if(!PyArg_ParseTuple(args, "Oi", &py_grid, &n_mas)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->n_mas= n_mas;

  Py_RETURN_NONE;
}

PyObject* py_grid_get_pk_normalisation(PyObject* self, PyObject* args)
{
  // _grid_get_pk_normalisation(_grid)
  // Return the degree of mass assighment
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("d", grid->pk_normalisation);
}

PyObject* py_grid_set_pk_normalisation(PyObject* self, PyObject* args)
{
  // _grid_set_nmas(_grid, pk_normalisation)
  PyObject *py_grid;
  double value;

  if(!PyArg_ParseTuple(args, "Od", &py_grid, &value)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->pk_normalisation= value;

  Py_RETURN_NONE;
}



PyObject* py_grid_fx_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  const int nd= 3;
  npy_intp nc= grid->nc;
  npy_intp ncz= 2*(nc/2 + 1);
  
  npy_intp dims[]= {nc, nc, nc};
  npy_intp strides[]= {(npy_intp) (sizeof(Float)*nc*ncz),
		       (npy_intp) (sizeof(Float)*ncz),
		       (npy_intp) (sizeof(Float))};
  return PyArray_New(&PyArray_Type, nd, dims, NPY_FLOAT_TYPE, strides,
  grid->fx, 0, NPY_ARRAY_WRITEABLE, 0);
}


PyObject* py_grid_fk_asarray(PyObject* self, PyObject* args)
{
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  const int nd= 3;
  const npy_intp nc= grid->nc;
  const int nckz= nc/2 + 1;
  npy_intp dims[]= {nc, nc, nckz};

  return PyArray_SimpleNewFromData(nd, dims, NPY_COMPLEX_TYPE, grid->fk);
}



PyObject* py_grid_clear(PyObject* self, PyObject* args)
{
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->clear();

  Py_RETURN_NONE;
}

PyObject* py_grid_compute_fluctuation(PyObject* self, PyObject* args)
{
  PyObject *py_grid_data, *py_grid_rand;
  if(!PyArg_ParseTuple(args, "OO", &py_grid_data, &py_grid_rand)) {
    return NULL;
  }

  Grid* const grid_data=
    (Grid*) PyCapsule_GetPointer(py_grid_data, "_Grid");
  py_assert_ptr(grid_data);

  Grid* const grid_rand=
    (Grid*) PyCapsule_GetPointer(py_grid_rand, "_Grid");
  py_assert_ptr(grid_rand);

  grid_compute_fluctuation(*grid_data, *grid_rand);

  Py_RETURN_NONE;
}


PyObject* py_grid_compute_fluctuation_homogeneous(PyObject* self,
						  PyObject* args)
{
  PyObject *py_grid_data;
  if(!PyArg_ParseTuple(args, "O", &py_grid_data)) {
    return NULL;
  }

  Grid* const grid_data=
    (Grid*) PyCapsule_GetPointer(py_grid_data, "_Grid");
  py_assert_ptr(grid_data);


  grid_compute_fluctuation_homogeneous(*grid_data);

  Py_RETURN_NONE;
}

PyObject* py_grid_load_fx_from_array(PyObject* self, PyObject* args)
{
  PyObject *py_grid, *py_array;

  if(!PyArg_ParseTuple(args, "OO", &py_grid, &py_array)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  Py_buffer a;
  if(PyObject_GetBuffer(py_array, &a, PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS | PyBUF_FULL_RO) == -1)
    return NULL;

  if(a.ndim != 3) {
    PyErr_SetString(PyExc_TypeError, "Expected a 3-dimensional array for grid");
    return NULL;
  }

  size_t nc= grid->nc;
  Py_ssize_t nc_check= static_cast<Py_ssize_t>(nc);
  if(a.shape[0] != nc_check || a.shape[1] != nc_check ||
     a.shape[2] != nc_check) {
    PyErr_SetString(PyExc_TypeError, "Expected a cubic grid.");
    return NULL;
  }

  if(strcmp(a.format, "d") != 0) {
    PyErr_SetString(PyExc_TypeError, "Expected an array of double.");
    return NULL;
  }
  
  npy_intp ncz= 2*(nc/2 + 1);

  double* p= (double*) a.buf;

  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      size_t index= (ix*nc + iy)*ncz;
      for(size_t iz=0; iz<nc; ++iz) {
	grid->fx[index++]= *p++;
      }
    }
  }
  Py_RETURN_NONE;
}


PyObject* py_grid_get_param(PyObject* self, PyObject* args)
{
  // _grid_get_param(_grid, param)
  // Return sum of weights, total_weight, w2_sum, nw2_sum
  PyObject* py_grid;
  char* param_name;

  if(!PyArg_ParseTuple(args, "Os", &py_grid, &param_name)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  string param(param_name);

  if(param == "shot_noise")
    return Py_BuildValue("d", grid->shot_noise);

  
  PyErr_SetString(PyExc_ValueError, "Unknown grid parameter");
  return NULL;
}

PyObject* py_grid_set_param_double(PyObject* self, PyObject* args)
{
  // _grid_set_sum(_grid, param, val)
  PyObject *py_grid;
  char* param_name;
  double val;

  if(!PyArg_ParseTuple(args, "Osd", &py_grid, &param_name, &val)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  string param(param_name);

  if(param == "shot_noise")
    grid->shot_noise= val;
  else {
    PyErr_SetString(PyExc_ValueError, "Unknown grid parameter");
    return NULL;
  }

  Py_RETURN_NONE;
}


PyObject* py_grid_resize_fourier(PyObject* self, PyObject* args)
{
  PyObject *py_grid, *py_grid_new;

  if(!PyArg_ParseTuple(args, "OO", &py_grid, &py_grid_new)) {
    return NULL;
  }

  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  Grid* const grid_new=
    (Grid*) PyCapsule_GetPointer(py_grid_new, "_Grid");
  py_assert_ptr(grid);

  grid_resize_fourier(grid, grid_new);

  Py_RETURN_NONE;
}

PyObject* py_grid_create_k(PyObject* self, PyObject* args)
{
  // Create a grid of k_i
  PyObject *py_grid;
  int axis;
  
  if(!PyArg_ParseTuple(args, "Oi", &py_grid, &axis)) {
    return NULL;
  }
  py_assert_ptr(0 <= axis && axis < 3);

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);
  
  const int nc= static_cast<int>(grid->nc);
  const Float boxsize= grid->boxsize;
  py_assert_ptr(nc > 0);
  py_assert_ptr(boxsize > 0.0);

  grid->mode= GridMode::fourier_space;

  const int nckz= nc/2 + 1;
  const Float fac= 2.0*M_PI/boxsize;
  complex<Float>* const d= grid->fk;

  Float k[3];
  for(int ix=0; ix<nc; ++ix) {
    k[0]= ix < nc/2 ? fac*ix : fac*(ix - nc);
    for(int iy=0; iy<nc; ++iy) {
      k[1]= iy < nc/2 ? fac*iy : fac*(iy - nc);
      for(int iz=0; iz<nckz; ++iz) {
	k[2]= fac*iz;
	
	size_t index= (ix + iy*nc)*nckz + iz;
	d[index] = k[axis];
      }
    }
  }
  Py_RETURN_NONE;
}

PyObject* py_grid_create_kmag(PyObject* self, PyObject* args)
{
  // Create a grid of |k|
  PyObject *py_grid;
  
  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  const int nc= static_cast<int>(grid->nc);
  const Float boxsize= grid->boxsize;
  py_assert_ptr(nc > 0);
  py_assert_ptr(boxsize > 0.0);

  grid->mode= GridMode::fourier_space;

  const int nckz= nc/2 + 1;
  const Float fac= 2.0*M_PI/boxsize;
  complex<Float>* d= grid->fk;

  for(int ix=0; ix<nc; ++ix) {
    Float kx= ix < nc/2 ? fac*ix : fac*(ix - nc);
    for(int iy=0; iy<nc; ++iy) {
      Float ky= iy < nc/2 ? fac*iy : fac*(iy - nc);
      for(int iz=0; iz<nckz; ++iz) {
	Float kz= fac*iz;
	
	size_t index= (ix + iy*static_cast<size_t>(nc))*nckz + iz;
	d[index] = sqrt(kx*kx + ky*ky + kz*kz);
      }
    }
  }
  Py_RETURN_NONE;
}

/*
PyObject* py_grid_set_k(PyObject* self, PyObject* args)
{
  // _grid_set_k(_grid)
  // Set k to the grid
  PyObject *py_grid;

  if(!PyArg_ParseTuple(args, "O", &py_grid)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->mode= GridMode::fourier_space;

  const Float fac= 2.0*M_PI/grid->boxsize;
  const int nc= static_cast<int>(grid->nc);
  const size_t nc_size= grid->nc;
  const size_t nckz= nc_size/2 + 1;
  
  const int iknq= nc/2;
  const int nkz= nc/2 + 1;

  complex<Float>* fk= grid->fk;
  
  for(int ix=0; ix<nc; ++ix) {
    int ikx= ix <= iknq ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      int iky= iy <= iknq ? iy : iy - nc;
      for(int iz=0; iz<nkz; ++iz) {
	int ikz= iz;
	size_t index= (ix*nc_size + iy)*nckz + iz;
	fk[index]= fac*sqrt(static_cast<double>(ikx*ikx + iky*iky + ikz*ikz));
      }
    }
  }

  Py_RETURN_NONE;
}
*/

PyObject* py_grid_set_mu2(PyObject* self, PyObject* args)
{
  // _grid_set_mu2(_grid, axis)
  // Set k to the grid
  PyObject *py_grid;
  int axis;
  
  if(!PyArg_ParseTuple(args, "Oi", &py_grid, &axis)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  grid->mode= GridMode::fourier_space;

  const int nc= static_cast<int>(grid->nc);
  const size_t nc_size= grid->nc;
  const size_t nckz= nc_size/2 + 1;
  
  const int iknq= nc/2;
  const int nkz= nc/2 + 1;

  complex<Float>* fk= grid->fk;
  int ik[3];
  
  for(int ix=0; ix<nc; ++ix) {
    ik[0]= ix <= iknq ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      ik[1]= iy <= iknq ? iy : iy - nc;

      for(int iz=0; iz<nkz; ++iz) {
	ik[2]= iz;
	Float k2= static_cast<Float>(ik[0]*ik[0] + ik[1]*ik[1] + ik[2]*ik[2]);
	Float mu2= static_cast<Float>(ik[axis]*ik[axis])/k2;
	size_t index= (ix*nc_size + iy)*nckz + iz;
	fk[index]= mu2;
      }
    }
  }

  fk[0]= 0.0;

  Py_RETURN_NONE;
}

PyObject* py_grid_set_power3d(PyObject* self, PyObject* args)
{
  // _grid_set_power3d(k, P, _grid)
  // Set P(k) to the grid
  PyObject *py_k, *py_P, *py_grid;
  
  if(!PyArg_ParseTuple(args, "OOO", &py_k, &py_P, &py_grid)) {
    return NULL;
  }

  vector<double> v_k, v_P;
  py_util_array_as_vector("k", py_k, v_k);
  py_util_array_as_vector("k", py_P, v_P, v_k.size());
  
  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  const double boxsize= grid->boxsize;
  const double fac= 2.0*M_PI/boxsize;

  // interpolate P(k)
  Interp interp(v_k, v_P);

  // Loop over all grid points
  const int nc= static_cast<int>(grid->nc);
  const size_t nc_size= grid->nc;
  const size_t nckz= nc_size/2 + 1;
  
  const int iknq= nc/2;
  const int nkz= nc/2 + 1;

  complex<Float>* fk= grid->fk;
  fk[0]= 0.0;

  int ik[3];
  for(int ix=0; ix<nc; ++ix) {
    ik[0]= ix <= iknq ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      ik[1]= iy <= iknq ? iy : iy - nc;

      int iz0= ix == 0 && iy == 0; // skip k = (0,0,0)
      for(int iz=iz0; iz<nkz; ++iz) {
	ik[2]= iz;
	Float k2= static_cast<Float>(ik[0]*ik[0] + ik[1]*ik[1] + ik[2]*ik[2]);
	Float k= fac*sqrt(k2);
	
	size_t index= (ix*nc_size + iy)*nckz + iz;
	fk[index]= interp(k);
      }
    }
  }

  grid->mode= GridMode::fourier_space;

  Py_RETURN_NONE;
}



PyObject* py_grid_write_binary_real(PyObject* self, PyObject* args)
{
  // _grid_write_binary_real(_grid, filename)
  // Write the grid to a binary file
  //
  // float_size: 4 write float
  // float_size: 8 write double
  // float_size: 0 write Float
  PyObject *py_grid, *bytes;
  char* filename;
  int float_size;
  Py_ssize_t len;

  if(!PyArg_ParseTuple(args, "OO&i", &py_grid, PyUnicode_FSConverter, &bytes,
		       &float_size)) {
    return NULL;
  }

  assert(float_size == 0 || float_size == 4 || float_size == 8);

  
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);
  py_assert_ptr(grid->mode == GridMode::real_space);
  
  FILE* fp= fopen(filename, "w");
  if(fp == 0) {
    char msg[256];
    sprintf(msg, "Unable to write binary file: %.200s", filename);
    PyErr_SetString(PyExc_FileNotFoundError, msg);
  }

  const size_t nc= grid->nc;
  const size_t ncz= 2*(nc/2+1);
  Float const* p= grid->fx;

  if(float_size == 0) {
    for(size_t i=0; i<nc*nc; ++i) {
      // write binary skipping the buffer region
      int ret= fwrite(p, sizeof(Float), nc, fp);
      assert(ret == static_cast<int>(nc));
      p += ncz;
    }
  }
  else if(float_size == 4)
    write_float_to_binary<float>(fp, p, nc);
  else if(float_size == 8)
    write_float_to_binary<double>(fp, p, nc);
  else
    assert(false);

  int ret= fclose(fp); py_assert_ptr(ret == 0);

  Py_DECREF(bytes);
  Py_RETURN_NONE;
}

PyObject* py_grid_write_vector_binary_real(PyObject* self, PyObject* args)
{
  // _grid_write_binary_vector_real(filename, _grid_x, _grid_y, grid_y)
  // Write the grid to a binary file
  //
  PyObject *bytes;
  char* filename;
  PyObject *py_grid_x, *py_grid_y, *py_grid_z;

  if(!PyArg_ParseTuple(args, "O&OOO", PyUnicode_FSConverter, &bytes,
		       &py_grid_x, &py_grid_y, &py_grid_z)) {
    return NULL;
  }

  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  Grid* const grid_x= (Grid*) PyCapsule_GetPointer(py_grid_x, "_Grid");
  Grid* const grid_y= (Grid*) PyCapsule_GetPointer(py_grid_y, "_Grid");
  Grid* const grid_z= (Grid*) PyCapsule_GetPointer(py_grid_z, "_Grid");
  py_assert_ptr(grid_x);
  py_assert_ptr(grid_y);
  py_assert_ptr(grid_z);
  py_assert_ptr(grid_x->mode == GridMode::real_space);
  py_assert_ptr(grid_y->mode == GridMode::real_space);
  py_assert_ptr(grid_z->mode == GridMode::real_space);
  
  FILE* fp= fopen(filename, "w");
  if(fp == 0) {
    char msg[256];
    sprintf(msg, "Unable to write binary file: %.200s", filename);
    PyErr_SetString(PyExc_FileNotFoundError, msg);
  }

  const size_t nc= grid_x->nc;
  const size_t ncz= 2*(nc/2+1);
  Float const* p_x= grid_x->fx;
  Float const* p_y= grid_y->fx;
  Float const* p_z= grid_z->fx;

  vector<float> v;
  v.reserve(3*nc);

  for(size_t i=0; i<nc*nc; ++i) {
    v.clear();
    for(size_t j=0; j<nc; ++j) {
      v.push_back(static_cast<float>(p_x[j]));
      v.push_back(static_cast<float>(p_y[j]));
      v.push_back(static_cast<float>(p_z[j]));
    }
    int ret= fwrite(v.data(), sizeof(float), 3*nc, fp);
    assert(ret == 3*static_cast<int>(nc));
    p_x += ncz;
    p_y += ncz;
    p_z += ncz;
  }  

  int ret= fclose(fp);
  py_assert_ptr(ret == 0);

  Py_DECREF(bytes);
  Py_RETURN_NONE;
}



PyObject* py_grid_copy(PyObject* self, PyObject* args)
{
  // _grid_copy(_grid)
  // Create a new grid object with same data
  PyObject *py_grid_src, *py_grid_dest;

  if(!PyArg_ParseTuple(args, "OO", &py_grid_src, &py_grid_dest)) {
    return NULL;
  }

  Grid const * const grid_src=
    (Grid const*) PyCapsule_GetPointer(py_grid_src, "_Grid");
  Grid * const grid_dest=
    (Grid*) PyCapsule_GetPointer(py_grid_dest, "_Grid");

  py_assert_ptr(grid_src);
  py_assert_ptr(grid_dest);

  grid_src->copy(grid_dest);

  Py_RETURN_NONE;
}

PyObject* py_grid_get_w2(PyObject* self, PyObject* args)
{
  PyObject *py_grid;
  int n;

  if(!PyArg_ParseTuple(args, "Oi", &py_grid, &n)) {
    return NULL;
  }

  py_assert_ptr(0 <= n && n <= 4);
  
  Grid const * const grid=
    (Grid const *) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  return Py_BuildValue("d", grid->w2_sum_n[n]);
}
