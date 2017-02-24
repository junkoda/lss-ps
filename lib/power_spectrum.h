#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H 1

#include <cmath>
#include "grid.h"

class PowerSpectrum {
 public:
  double kmin, kmax, dk;
  int* nmode_hist;
  double *k_hist, *P0_hist, *P2_hist, *P4_hist;
  int n;

  PowerSpectrum(const double kmin_, const double kmax_, const double dk_);
  ~PowerSpectrum();

  void clear();
  void normalise();
  

  void add(const double k, const double P0,
	   const double P2, const double P4) {
    int i= (int) floor((k - kmin)*dk_inv);
    if(0 <= i && i < n) {
      nmode_hist[i]++;
      k_hist[i] += k;
      P0_hist[i] += P0;
      P2_hist[i] += P2;
      P4_hist[i] += P4;
    }
  }
 private:
  double dk_inv;
};

void power_spectrum_compute_multipoles(Grid const * const grid,
				       PowerSpectrum* const ps,
				       const bool subtract_shotnoise);

/*
void power_spectrum_compute_multipoles_raw(Grid const * const grid,
					   PowerSpectrum* const ps);

void power_spectrum_compute_multipoles(Grid const * const grid,
				       PowerSpectrum* const ps,
				       const bool subtract_shotnoise,
				       const Float neff);

void power_spectrum_compute_multipoles_interlacing2(
			    Grid const * const grid,
			    PowerSpectrum* const ps,
			    const bool subtract_shotnoise);

void power_spectrum_compute_shotnoise(Grid const * const grid,
				      PowerSpectrum* const ps);
*/

#endif
