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
#include "py_yamamoto.h"
#include "py_performance.h"
#include "py_mean_density.h"
#include "py_mean_density_adaptive.h"
#include "py_etc.h"
#include "py_gadget_file.h"
#include "py_cola_binary.h"
#include "py_model.h"
#include "py_rr.h"
#include "py_lognormal.h"
#include "py_window.h"

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

  {"_cosmology_growth_D", py_cosmology_growth_D, METH_VARARGS, "growth factor D"},
  {"_cosmology_growth_f", py_cosmology_growth_f, METH_VARARGS, "growth rate f"},
   
  {"_grid_alloc", py_grid_alloc, METH_VARARGS,
   "_grid_alloc(nc) allocate a new grid object"},
  {"_grid_fft", py_grid_fft, METH_VARARGS,
   "_grid_fft(_grid); apply FFT to the grid"},
  {"_grid_fft_inverse", py_grid_fft_inverse, METH_VARARGS,
   "_grid_fft_inverser(_grid); apply inverse FFT to the grid"},
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
  {"_grid_get_pk_normalisation", py_grid_get_pk_normalisation, METH_VARARGS,
   "set pk_normalisation"},
  {"_grid_set_pk_normalisation", py_grid_set_pk_normalisation, METH_VARARGS,
   "set pk_normalisation"},
  {"_grid_get_param", py_grid_get_param, METH_VARARGS,
   "get a grid parameter"},
  {"_grid_set_param_double", py_grid_set_param_double, METH_VARARGS,
   "set a floating point parameter"},
  {"_grid_resize_fourier", py_grid_resize_fourier, METH_VARARGS,
   "_grid_resize_fourier(_grid, _grid_new)"},

  {"_grid_fx_asarray", py_grid_fx_asarray, METH_VARARGS,
   "return real-space grid as an np.array"},
  {"_grid_fk_asarray", py_grid_fk_asarray, METH_VARARGS,
   "return Fourier-space grid as an np.array"},
  {"_grid_load_fx_from_array", py_grid_load_fx_from_array, METH_VARARGS,
   "fill grid data from array"},
  
  {"_grid_clear", py_grid_clear, METH_VARARGS,
   "_grid_clear(_grid); clear the array with 0"},
  {"_grid_copy", py_grid_copy,  METH_VARARGS,
   "_grid_copy(_grid_src, _grid_dest); copy src grid to dest grid"},   
  {"_grid_compute_fluctuation", py_grid_compute_fluctuation, METH_VARARGS,
   "_grid_compute_fluctuation(_grid_data, _grid_rand)"},
  {"_grid_compute_fluctuation_homogeneous",
   py_grid_compute_fluctuation_homogeneous, METH_VARARGS,
   "_grid_compute_fluctuation_homogeneous(_grid_data)"},

  {"_grid_create_k", py_grid_create_k, METH_VARARGS,
   "_grid_create_k(_grid, nc, boxsize, axis)"},
  {"_grid_create_kmag", py_grid_create_kmag, METH_VARARGS,
   "_grid_create_kmag(_grid, nc, boxsize)"},
  {"_grid_set_mu2", py_grid_set_mu2, METH_VARARGS,
   "_grid_set_mu2(_grid, axis); set mu2 = (k[axis]/|k|)^2"},
  {"_grid_set_power3d", py_grid_set_power3d, METH_VARARGS,
   "_grid_set_power3d(k, P, _grid)"},
  
  {"_grid_write_binary_real", py_grid_write_binary_real, METH_VARARGS,
   "_grid_write_binary_real(_grid, filename, float_size)"},
  {"_grid_write_vector_binary_real",
   py_grid_write_vector_binary_real, METH_VARARGS,
   "_grid_write_vector_binary_real(filename, _grid_x, _grid_y, grid_z)"},
  
  {"_mass_assignment_from_array", py_mass_assignment_from_array, METH_VARARGS,
   "assign mass to a grid from an array"},
  {"_mass_assignment_variations", py_mass_assignment_variations, METH_VARARGS,
   "assign mass to a grid from an array (serial or atomic)"},
  {"_mass_assignment_correct_mas", py_mass_assignment_correct_mas, METH_VARARGS,
   "correct mass assignment window function"},

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
  {"_power_spectrum_Pl_asarray", py_power_spectrum_Pl_asarray, METH_VARARGS,
   "_power_spectrum_Pl_asarray(_ps, l) P[l] as array"},
  {"_power_spectrum_compute_plane_parallel",
   py_power_spectrum_compute_plane_parallel, METH_VARARGS,
   "_power_spectrum_compute_plane_parallel(k_min, k_max, dk, _grid_delta, subtract_shotnoise, correct_mas)"},
  {"_power_spectrum_compute_yamamoto",
   py_power_spectrum_compute_yamamoto, METH_VARARGS,
   "_power_spectrum_compute_yamamoto(k_min, k_max, dk, grid, grid2, grid4, subtract_shotnoise, correct_mas)"},
  {"_power_spectrum_compute_yamamoto_odd",
   py_power_spectrum_compute_yamamoto_odd, METH_VARARGS,
   "_power_spectrum_compute_yamamaoto_odd(k_min, k_max, dk, _grid, _grid1, _grid3, subtract_shotnoise, correct_mas)"},


  {"_power_spectrum_shotnoise",py_power_spectrum_shotnoise, METH_VARARGS,
   "_power_spectrum_shotnoise(_ps)"},
  {"_power_spectrum_compute_power_multipoles",
   py_power_spectrum_compute_power_multipoles, METH_VARARGS,
   "compute multipoles of 3D power spectrum grid"},
  {"_power_spectrum_compute_discrete_multipoles",
   py_power_spectrum_compute_discrete_multipoles, METH_VARARGS,
   "compute discrete legendre multipoles of delta(k) or P(k)"},

  {"_power_spectrum_compute_2d_average", py_power_spectrum_compute_2d_average,
   METH_VARARGS, "compute 2d average of a give 3D grid"},
  
  {"_yamamoto_compute_moment_x", py_yamamoto_compute_moment_x, METH_VARARGS,
   "_yamamoto_compute_moment_x"},  
  {"_yamamoto_compute_moment_k", py_yamamoto_compute_moment_k, METH_VARARGS,
   "_yamamoto_compute_moment_k"},   

  {"_performance_xyz_sum", py_performance_xyz_sum, METH_VARARGS,
   "compute sum of an xyz array"},

  {"_mean_density_from_grid", py_mean_density_from_grid, METH_VARARGS,
   "_mean_density_from_grid(_grid, xyz, nbar)"},

  {"_etc_gaussian_smoothing", py_etc_gaussian_smoothing, METH_VARARGS,
   "_etc_gaussian_smoothing(_grid, sigma_smoothing)"},

  {"_kdpoints_alloc", py_kdpoints_alloc, METH_VARARGS,
   "_kdpoints_alloc()"},
  {"_kdpoints_from_array", py_kdpoints_from_array, METH_VARARGS,
   "_kdpoints_from_array(_kdpoints, xyz, w)"},
  {"_kdpoints_density_as_array", py_kdpoints_density_as_array, METH_VARARGS,
   "_kdpoints_density_as_array(_points, nbar)"},
  {"_kdpoints_len", py_kdpoints_len, METH_VARARGS,
   "_kdpoints_len(_points)"},

  {"_kdtree_alloc", py_kdtree_alloc, METH_VARARGS,
   "_kdtree_alloc()"},
  {"_kdtree_compute_rk", py_kdtree_compute_rk, METH_VARARGS,
   "_kdtree_compute_rk(_kdtree, knbr)"},
  {"_kdtree_update_node_statistics", py_kdtree_update_node_statistics,
   METH_VARARGS, "_kdtree_update_node_statistics(_kdtree)"},
  {"_mean_density_adaptive_estimate", py_mean_density_adaptive_estimate,
   METH_VARARGS, "_mean_density_adaptive_estimate(_kdtree, _points, knbr)"},
  {"_mean_density_average_estimate", py_mean_density_average_estimate,
   METH_VARARGS, "_mean_density_average_estimate(_kdtree, xyz, knbr, nbar)"},   

  {"_model_linear_alloc", py_model_linear_alloc, METH_VARARGS,
   "_model_linear_alloc(filename, omega_m, z, b, redshift_space, sigma_v)"},
  {"_model_call",py_model_call, METH_VARARGS,
   "_model_call(_model, k, mu)"},
  {"_model_compute_discrete_multipoles", py_model_compute_discrete_multipoles,
   METH_VARARGS, "_model_compute_descrete_multipoles"},
  {"_model_apply_window_3d", py_model_apply_window_3d,
   METH_VARARGS, "_model_apply_window_3d(_model, _grid, k)"},
  {"_model_create_grid", py_model_create_grid, METH_VARARGS,
   "_model_create_grid"},

  {"_gadget_file_read_np", py_gadget_file_read_np, METH_VARARGS,
   "_gadget_file_read_np(filename)"},
  {"_gadget_file_read_header", py_gadget_file_read_header, METH_VARARGS,
   "_gadget_file_read_header"},
  {"_gadget_file_alloc", py_gadget_file_alloc, METH_VARARGS,
   "_gadget_file_alloc(filename)"},
  {"_gadget_file_open", py_gadget_file_open, METH_VARARGS,
   "_gadget_file_open(_gf)"},
  {"_gadget_file_close", py_gadget_file_close, METH_VARARGS,
   "_gadget_file_close(_gf)"},
  {"_gadget_file_read", py_gadget_file_read, METH_VARARGS,
   "_gadget_file_read(_gf, components, ibegin, iend, a)"},

  {"_cola_binary_load_particles_header",
   py_cola_binary_load_particles_header, METH_VARARGS,
   "_cola_binary_load_particles_header(filename)"},
  {"_cola_binary_load_particles",
   py_cola_binary_load_particles, METH_VARARGS,
   "_cola_binary_load_particles(filename, a)"},

  {"_rr_compute_multipoles", py_rr_compute_multipoles, METH_VARARGS,
   "compute RR multipoles for window function convolution"},
  {"_rr_asarray", py_rr_asarray, METH_VARARGS,
   "_rr_asarray(_rr, l); get rr[l] as an array"},

  {"_lognormal_convert_P_to_delta_k",
   py_lognormal_convert_P_to_delta_k, METH_VARARGS,
   "_logrnomal_convert_P_to_delta_k(_grid, seed, fix_amplitude"},

  {"_window_cquag_bessel_transform", 
   py_window_cquag_bessel_transform, METH_VARARGS,
   "_window_cquag_bessel_transform()"},

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
  py_rr_module_init();
  //py_gadget_file_module_init();
  //py_mean_density_adaptive_module_init();
  
  return PyModule_Create(&module);
}
