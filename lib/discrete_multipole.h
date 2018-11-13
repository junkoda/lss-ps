#ifndef DISCRETE_MULTIPOLE_H
#define DISCRETE_MULTIPOLE_H 1

#include <vector>

void discrete_multipole_compute_legendre(const double k_min, const double k_max, const double dk, const Float boxsize, std::vector<Float>& coef);

#endif
