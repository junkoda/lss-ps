#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <cstdlib>
void cosmology_init(const double omega_m_, const double z_max_, const size_t n);
void cosmology_free();

double cosmology_omega_m();
double cosmology_compute_comoving_distance(const double a);
double cosmology_distance_redshift(const double d);
double cosmology_redshift_distance(const double z);




#endif
