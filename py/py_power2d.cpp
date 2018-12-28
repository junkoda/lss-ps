#include "grid.h"
#include "power2d.h"
#include "py_util.h"
#include "py_buffer.h"
#include "py_power2d.h"

PyObject* py_power2d_compute(PyObject* self, PyObject* args)
{
  // _power2d_compute(grid1, grid2, real_imag,
  //                  k_min, dk, nk, nmu,
  //                  shot_noise, correct_mas,
  //                  nmodes_array, k_array, mu_array, ps2d_array) 

  PyObject *py_grid1, *py_grid2;
  int real_imag;
  double k_min, dk;
  int nk, nmu;
  double shot_noise;
  int correct_mas, los;
  PyObject *py_nmodes, *py_k, *py_mu, *py_ps2d;
  
  if(!PyArg_ParseTuple(args, "OOiddiidiiOOOO",
		       &py_grid1, &py_grid2, &real_imag,
		       &k_min, &dk, &nk, &nmu,
		       &shot_noise, &correct_mas, &los,
		       &py_nmodes, &py_k, &py_mu, &py_ps2d))
    return NULL;

  Grid* const grid1= (Grid*) PyCapsule_GetPointer(py_grid1, "_Grid");
  py_assert_ptr(grid1);

  Grid* const grid2= (Grid*) PyCapsule_GetPointer(py_grid2, "_Grid");
  py_assert_ptr(grid2);


  //
  // Decode array information
  //
  //Py_buffer buf_nmodes, buf_k, buf_mu, buf_ps2d;
  Buffer buf_nmodes("nmodes2d", py_nmodes);
  Buffer buf_k("k2d", py_k);
  Buffer buf_mu("mu2d", py_mu);
  Buffer buf_ps2d("ps2d", py_ps2d);

  try {
    Py_ssize_t len= buf_nmodes.shape(0);
    
    buf_nmodes.assert_shape(1, len);
    buf_k.assert_shape(1, len);
    buf_mu.assert_shape(1, len);
    buf_ps2d.assert_shape(1, len);

    py_assert_ptr(nk*nmu == static_cast<int>(len));
  }
  catch(...) {
    return NULL;
  }

  try {
    power2d_compute(grid1, grid2,real_imag,
		    k_min, dk, nk, nmu,
		    shot_noise, correct_mas, los,
		    buf_nmodes.buf<double>(),
		    buf_k.buf<double>(),
		    buf_mu.buf<double>(),
		    buf_ps2d.buf<double>());
  }
  catch(...) {
    return NULL;
  }

  /*
  PyBuffer_Release(&buf_nmodes);
  PyBuffer_Release(&buf_k);
  PyBuffer_Release(&buf_mu);
  PyBuffer_Release(&buf_ps2d);
  */

  Py_RETURN_NONE;
}
