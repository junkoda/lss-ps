#ifndef MULTIPOLE_H
#define MULTIPOLE_H 1

#include "power_spectrum.h"
#include "grid.h"

PowerSpectrum*
multipole_compute_monopole(const double k_min, const double k_max,
			   const double dk, const int nmu,
			   Grid const * const grid,
			   const bool subtract_shotnoise,
			   const bool correct_mas);

PowerSpectrum*
multipole_compute_plane_parallel(const double k_min, const double k_max,
				 const double dk, const int nmu,
				 Grid const * const grid,
				 const bool subtract_shotnoise,
				 const bool correct_mas,
				 const int line_of_sight);

PowerSpectrum*
multipole_compute_yamamoto_scoccimarro(const double k_min, const double k_max,
			   const double dk,
			   Grid const * const grid,
			   Grid const * const grid2,
			   const bool subtract_shotnoise,
			   const bool correct_mas);

PowerSpectrum*
multipole_compute_yamamoto_bianchi(const double k_min, const double k_max,
			   const double dk,
			   Grid const * const grid,
			   Grid const * const grid2,
			   Grid const * const grid4,
			   const bool subtract_shotnoise,
			   const bool correct_mas);



#endif
