#include <vector>
#include "mean_density_adaptive.h"
#include "py_mean_density_adaptive.h"
#include "py_assert.h"

/*
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#ifdef DOUBLEPRECISION
#define NPY_FLOAT_TYPE NPY_DOUBLE
#else
#define NPY_FLOAT_TYPE NPY_FLOAT
#endif


PyMODINIT_FUNC
py_mean_density_adaptive_module_init()
{
  import_array();

  return NULL;
}
*/

using namespace std;

// Create a vector of KDPoint from pointers of float_type (double or float)
template <typename float_type>
void construct_vector(float_type const * xyz,
		      const size_t xyz_stride,
		      float_type const * weight,
		      const size_t weight_stride,
		      const size_t np,
		      vector<KDPoint>& v)
{
  v.clear();
  v.reserve(np);

  KDPoint p;
  p.w = 1;
  p.n_local = 0;
  p.n_average = 0;
  
  for(size_t i=0; i<np; ++i) {
    p.idx= i; // original index
    p.x[0] = xyz[0];
    p.x[1] = xyz[1];
    p.x[2] = xyz[2];
    xyz= (float_type*) ((char*) xyz + xyz_stride);
    
    if(weight) {
      p.w = *weight;
      weight = (float_type*) ((char*) weight + weight_stride);
    }
    v.push_back(p);
  }
}

static void py_kdtree_free(PyObject *obj);
static void py_kdpoints_free(PyObject *obj);

//
// KDPoints
//   wrapper for vector<KDPoint>
//

PyObject* py_kdpoints_alloc(PyObject* self, PyObject* args)
{
  // _kdpoints_alloc_empty()

  // Create a new empty vector<KDPoint>
  vector<KDPoint>* const pv= new vector<KDPoint>();

  return PyCapsule_New(pv, "_KDPoints", py_kdpoints_free);
}

void py_kdpoints_free(PyObject *obj)
{
  // Delete the _kdpoints object
  // Called automatically by Python
   vector<KDPoint>* const pv=
    (vector<KDPoint>*) PyCapsule_GetPointer(obj, "_KDPoints");
  py_assert_void(pv);

  delete pv;
}


PyObject* py_kdpoints_from_array(PyObject* self, PyObject* args)
{
  // Construct vector<KDPoint> from arrays of xyz and weights
  // _kdpoints_alloc_from_array(_kdpoints, xyz, weight)
  PyObject *py_kdpoints, *py_xyz, *py_weight;

  if(!PyArg_ParseTuple(args, "OOO", &py_kdpoints, &py_xyz, &py_weight)) {
    return NULL;
  }

  vector<KDPoint>* const pv= 
    (vector<KDPoint>*) PyCapsule_GetPointer(py_kdpoints, "_KDPoints");
  py_assert_ptr(pv);

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

  // weight array
  Py_buffer weight;
  Py_ssize_t zero[]= {(Py_ssize_t) 0};
  if(py_weight == Py_None) {
    weight.buf= nullptr;
    weight.suboffsets= nullptr;
    weight.strides= zero;
  }
  else {
    if(PyObject_GetBuffer(py_weight, &weight, PyBUF_FORMAT | PyBUF_FULL_RO)
       == -1) return NULL;

    if(weight.ndim != 1) {
      PyErr_SetString(PyExc_TypeError,
		      "Expected a 1-dimensional array for weight");
      return NULL;
    }

    if(weight.shape[0] != xyz.shape[0]) {
      PyErr_SetString(PyExc_TypeError, "Length of xyz and weight arrays differ");
      return NULL;
    }

    if(strcmp(xyz.format, weight.format) != 0) {
      PyErr_SetString(PyExc_TypeError,
	 "data type (double/single precision) of xyz and weight are not the same");
      return NULL;
    }
  }

  if(xyz.suboffsets || weight.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "_mean_density_adaptive cannot handle array with suboffsets");
    return NULL;
  }
    
  if(strcmp(xyz.format, "d") == 0) {
    if(xyz.strides[1] != sizeof(double)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    construct_vector<double>((double*) xyz.buf,
			     xyz.strides[0],
			     (double*) weight.buf,
			     weight.strides[0],
			     xyz.shape[0],
			     *pv);
  }
  else if(strcmp(xyz.format, "f") == 0) {
    if(xyz.strides[1] != sizeof(float)) {
      PyErr_SetString(PyExc_TypeError,
		    "xyz array is expected to be contigeous in 2nd direction");
      return NULL;
    }

    construct_vector<float>((float*) xyz.buf,
			    xyz.strides[0],
			    (float*) weight.buf,
			    weight.strides[0],
			    xyz.shape[0],
			    *pv);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "Expected an array of floats or doubles");
    return NULL;
  }

  Py_RETURN_NONE;
}

PyObject* py_kdpoints_len(PyObject* self, PyObject* args)
{
  // _kdpoints_len(_points)
  PyObject *py_kdpoints;

  if(!PyArg_ParseTuple(args, "O", &py_kdpoints)) {
    return NULL;
  }

  vector<KDPoint>* const pv= 
    (vector<KDPoint>*) PyCapsule_GetPointer(py_kdpoints, "_KDPoints");
  py_assert_ptr(pv);

  return Py_BuildValue("K", (unsigned long long) pv->size());
}
  
PyObject* py_kdpoints_density_as_array(PyObject* self, PyObject* args)
{
  // Copy points.density to the given array
  // _kdpoints_alloc_from_array(_kdpoints, nbar)
  PyObject *py_kdpoints, *py_nbar;

  if(!PyArg_ParseTuple(args, "OO", &py_kdpoints, &py_nbar)) {
    return NULL;
  }

  vector<KDPoint> const * const pv= 
    (vector<KDPoint> const *) PyCapsule_GetPointer(py_kdpoints, "_KDPoints");
  py_assert_ptr(pv);
  const vector<KDPoint>& v= *pv;

  //
  // Decode array information
  //

  // output nbar array
  Py_buffer nbar;
  if(PyObject_GetBuffer(py_nbar, &nbar, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;
  
  if(nbar.ndim != 1) {
    PyErr_SetString(PyExc_TypeError, "Expected a 1-dimensional array for nbar");
    return NULL;
  }

  if(nbar.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
	    "_mean_density_adaptive cannot handle array with suboffsets");
    return NULL;
  }

  if(static_cast<size_t>(nbar.shape[0]) != pv->size()) {
    PyErr_SetString(PyExc_TypeError, "Length of points and nbar array differ");
    return NULL;
  }

  if(strcmp(nbar.format, "d") != 0) {
    PyErr_SetString(PyExc_TypeError, "Expeced an array of double for nbar");
    return NULL;
  }

  const size_t n= pv->size();
  double* const out= (double*) nbar.buf;
  
  for(size_t i=0; i<n; ++i) {
    // kdtree construction modifies the order of particles
    // returns density in the original order
    index_t idx= v[i].idx;
    assert(0 <= idx && idx < static_cast<index_t>(n));
    out[idx]= v[i].n_local;
  }

  cerr << "py 2852 " << v[2852].n_local << " " << v[2852].idx << endl;

  Py_RETURN_NONE;
}


//
// KDTree
//
PyObject* py_kdtree_alloc(PyObject* self, PyObject* args)
{
  // _kdtree_alloc(_kdpoints, knbr)
  //     xyz:    array of xyz
  //     weight: array of weights, can be None

  // Create a new KDTree object

  PyObject *py_points;
  int quota;

  if(!PyArg_ParseTuple(args, "Oi", &py_points, &quota)) {
    return NULL;
  }

  vector<KDPoint>* const pv=
    (vector<KDPoint>*) PyCapsule_GetPointer(py_points, "_KDPoints");
  py_assert_ptr(pv);

  KDTree* const kdtree= new KDTree(*pv, quota);

  return PyCapsule_New(kdtree, "_KDTree", py_kdtree_free);
}

void py_kdtree_free(PyObject *obj)
{
  // Delete the KDTree object
  // Called automatically by Python
  KDTree* const kdtree=
    (KDTree*) PyCapsule_GetPointer(obj, "_KDTree");
  py_assert_void(kdtree);

  delete kdtree;
}

//
// Density estimation
//
PyObject* py_mean_density_adaptive_estimate(PyObject* self, PyObject* args)
{
  // _mean_density_adaptive_estimate(_kdtree, _points, knbr)
  
  // Args:
  //     _kdtree (_KDTree)
  //     _points (_KDPoints)
  //     knbr (int): number of neighbors for the density estimation
  // Output:
  //     points.density

  PyObject *py_kdtree, *py_points;
  int knbr;

  if(!PyArg_ParseTuple(args, "OOi", &py_kdtree, &py_points, &knbr))
    return NULL;

  KDTree* const kdtree=
    (KDTree*) PyCapsule_GetPointer(py_kdtree, "_KDTree");
  py_assert_ptr(kdtree);

  vector<KDPoint>* const pv= 
    (vector<KDPoint>*) PyCapsule_GetPointer(py_points, "_KDPoints");
  py_assert_ptr(pv);

  kdtree->adaptive_kernel_density(*pv, knbr);

  Py_RETURN_NONE;
}

