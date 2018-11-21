#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H 1

#include <iostream>
#include <valarray>
#include <cmath>
#include "grid.h"

class PowerSpectrum {
 public:
  PowerSpectrum();
  PowerSpectrum(const int nk);
  
  PowerSpectrum& operator+=(const PowerSpectrum& ps);

  int n;
  double shot_noise;
  
  std::valarray<double> nmodes, nmodes2d;
  std::valarray<double> k, p0, p1, p2, p3, p4;
  //std::valarray<double> p2d;
 private:
  PowerSpectrum(const PowerSpectrum&);
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
