#ifndef POWER2D_H
#define POWER2D_H 1

#include "grid.h"

void power2d_compute(Grid const * const grid1,
		     Grid const * const grid2,
		     const int real_imag,
		     const double k_min,
		     const double dk,
		     const int nk, const int nmu,
		     const double shot_noise,
		     const bool correct_mas,
		     const int los,
		     double * const nmodes_out,
		     double * const k_out,
		     double * const mu_out,
		     double * const P_out);

#endif
