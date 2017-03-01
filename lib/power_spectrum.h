#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H 1

#include <iostream>
#include <cmath>
#include "grid.h"

class PowerSpectrum {
 public:
  // k_min, k_max, dk in Mpc/h
  // ik_min, ik_max are in units of k_fundamental = 2pi/boxsize
  const double kmin, dk;
  double kmax;
  const double k_fundamental;
  int* nmode_hist;
  double *k_hist, *P0_hist, *P2_hist, *P4_hist;
  const int n;

  PowerSpectrum(const double kmin_, const double kmax_, const double dk_,
		const double boxsize);
  ~PowerSpectrum();

  void clear();
  void normalise();

  void add(const double ik, const double P0,
	   const double P2, const double P4) {
    // k is in units of fundamental frequency
    //int i= (int) floor((ik - ik_min)*idk_inv);
    int i= (int) floor(ik*idk_inv); // DEBUG!!!

    //kstd::cerr << ik << " " << i << " " << idk_inv << std::endl;
    
    if(0 <= i && i < n) {
      nmode_hist[i]++;
      k_hist[i] += ik;
      P0_hist[i] += P0;
      P2_hist[i] += P2;
      P4_hist[i] += P4;
    }
  }

  int n_modes(const int i) {
    return nmode_hist[i];
  }
  double k(const int i) {
    // DEBUG
    return k_fundamental*k_hist[i];
  }
  double P0(const int i) {
    return P0_hist[i]; ///nmode_hist[i];
  }
  double P2(const int i) {
    return P2_hist[i];  ///nmode_hist[i];
  }
  double P4(const int i) {
    return P4_hist[i]; ///nmode_hist[i];
  }

 private:
  const double ik_min, idk_inv;
  double ik_max;
  //double dik_inv;
};

void power_spectrum_compute_multipoles(Grid const * const grid,
				       PowerSpectrum* const ps,
				       const bool subtract_shotnoise,
				       const bool correct_mas);

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
