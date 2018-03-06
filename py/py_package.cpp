//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_config.h"
#include "py_msg.h"
#include "py_cosmology.h"
#include "py_grid.h"
#include "py_mass_assignment.h"
#include "py_interlacing.h"
#include "py_power_spectrum.h"
#include "py_performance.h"
#include "py_mean_density.h"

using namespace std;

static PyMethodDef methods[] = {
  {"_config_sizeof_float", py_config_sizeof_float, METH_VARARGS,
   "return sizeof(Float)"},

  {"_msg_get_loglevel", py_msg_get_loglevel, METH_VARARGS,
   "return current loglevel"},
  {"_msg_set_loglevel", py_msg_set_loglevel, METH_VARARGS,
   "_msg_set_loglevel(lv); set loglevel"},

  {"_cosmology_init", py_cosmology_init, METH_VARARGS,
   "_cosmology_init(omega_m, z_max, n)"},
  {"_cosmology_distance_redshift", py_cosmology_distance_redshift,
   METH_VARARGS, "_cosmology_distance_redshift(d_array)"},
  {"_cosmology_compute_comoving_distance",
   py_cosmology_compute_comoving_distance, METH_VARARGS,
   "_cosmology_compute_comoving_distance"},
   
  {"_grid_alloc", py_grid_alloc, METH_VARARGS,
   "_grid_alloc(nc) allocate a new grid object"},
  {"_grid_fft", py_grid_fft, METH_VARARGS,
   "_grid_fft(_grid); apply FFT to the grid"},
  {"_grid_nc", py_grid_nc, METH_VARARGS, "return nc"},
  {"_grid_get_mode", py_grid_get_mode, METH_VARARGS, "return grid mode"},
  {"_grid_set_mode", py_grid_set_mode, METH_VARARGS, "set grid mode"},
  {"_grid_get_boxsize", py_grid_get_boxsize, METH_VARARGS, "return boxsize"},
  {"_grid_set_boxsize", py_grid_set_boxsize, METH_VARARGS, "set boxsize"},
  {"_grid_get_x0", py_grid_get_x0, METH_VARARGS, "return x0"},
  {"_grid_set_x0", py_grid_set_x0, METH_VARARGS, "set x0"},
  {"_grid_get_offset", py_grid_get_offset, METH_VARARGS, "return offset"},
  {"_grid_set_offset", py_grid_set_offset, METH_VARARGS, "set offset"},
  {"_grid_get_sums", py_grid_get_sums, METH_VARARGS, "get sums"},
  {"_grid_set_sums", py_grid_set_sums, METH_VARARGS, "set sums"},
  {"_grid_get_nmas", py_grid_get_nmas, METH_VARARGS, "get n_mas"},
  {"_grid_set_nmas", py_grid_set_nmas, METH_VARARGS, "set n_mas"},

  {"_grid_fx_asarray", py_grid_fx_asarray, METH_VARARGS,
   "return real-space grid as an np.array"},
  {"_grid_fk_asarray", py_grid_fk_asarray, METH_VARARGS,
   "return Fourier-space grid as an np.array"},
  {"_grid_load_fx_from_array", py_grid_load_fx_from_array, METH_VARARGS,
   "fill grid data from array"},
  
  {"_grid_clear", py_grid_clear, METH_VARARGS,
   "_grid_clear(_grid); clear the array with 0"},
  {"_grid_compute_fluctuation", py_grid_compute_fluctuation, METH_VARARGS,
   "_grid_compute_fluctuation(_grid_data, _grid_rand)"},
  {"_grid_compute_fluctuation_homogeneous",
   py_grid_compute_fluctuation_homogeneous, METH_VARARGS,
   "_grid_compute_fluctuation_homogeneous(_grid_data)"},

  {"_mass_assignment_from_array", py_mass_assignment_from_array, METH_VARARGS,
   "assign mass to a grid from an array"},
  {"_mass_assignment_variations", py_mass_assignment_variations, METH_VARARGS,
   "assign mass to a grid from an array (serial or atomic)"},

  {"_interlacing", py_interlacing, METH_VARARGS,
   "_interlacing(_grid1, _grid2)"},

  {"_power_spectrum_len",py_power_spectrum_len, METH_VARARGS,
   "_power_spectrum_len(_ps)"},
  {"_power_spectrum_k_asarray", py_power_spectrum_k_asarray, METH_VARARGS,
   "k array"},
  {"_power_spectrum_nmodes_asarray", py_power_spectrum_nmodes_asarray,
   METH_VARARGS, "number of modes array"},
  {"_power_spectrum_P0_asarray", py_power_spectrum_P0_asarray, METH_VARARGS,
   "P0 array"},   
  {"_power_spectrum_P2_asarray", py_power_spectrum_P2_asarray, METH_VARARGS,
   "P2 array"},   
  {"_power_spectrum_P4_asarray", py_power_spectrum_P4_asarray, METH_VARARGS,
   "P4 array"},
  {"_power_spectrum_compute_plane_parallel",
   py_power_spectrum_compute_plane_parallel, METH_VARARGS,
   "_power_spectrum_compute_plane_parallel(k_min, k_max, dk, nmu, _grid_delta, subtract_shotnoise, correct_mas)"},

  {"_performance_xyz_sum", py_performance_xyz_sum, METH_VARARGS,
   "compute sum of an xyz array"},

  {"_mean_density_from_grid", py_mean_density_from_grid, METH_VARARGS,
   "_mean_density_from_grid(_grid, xyz, nbar)"},
  
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_lssps", // name of this module
  "A package for power spectrum cosmology:large-scale structure", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__lssps(void) {
  py_grid_module_init();
  py_power_spectrum_module_init();    
  
  return PyModule_Create(&module);
}
