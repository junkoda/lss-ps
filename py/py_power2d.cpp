#include "grid.h"
#include "power2d.h"
#include "py_util.h"
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
  Py_ssize_t zero[]= {(Py_ssize_t) 0};
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
  Py_buffer buf_nmodes, buf_k, buf_mu, buf_ps2d;

  py_util_decode_array("nmodes2d", py_nmodes, &buf_nmodes, 0, 0, false);
  Py_ssize_t len= buf_nmodes.shape[0];
  py_util_decode_array("k2d", py_k, &buf_k, len, 0, false);
  py_util_decode_array("mu2d", py_mu, &buf_mu, len, 0, false);
  py_util_decode_array("ps2d", py_ps2d, &buf_ps2d, len, 0, false);

  py_assert_ptr(nk*nmu == static_cast<int>(len));

  power2d_compute(grid1, grid2,real_imag,
		  k_min, dk, nk, nmu,
		  shot_noise, correct_mas, los,
		  (double *) buf_nmodes.buf,
		  (double *) buf_k.buf,
		  (double *) buf_mu.buf,
		  (double *) buf_ps2d.buf);

  PyBuffer_Release(&buf_nmodes);
  PyBuffer_Release(&buf_k);
  PyBuffer_Release(&buf_mu);
  PyBuffer_Release(&buf_ps2d);

  Py_RETURN_NONE;
}
