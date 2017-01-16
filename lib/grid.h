#ifndef GRID_H
#define GRID_H 1

#include <fftw3.h>
#include "config.h"

// FFTW() macro adds fftw_ or fftwf_ prefix depending on DOUBLEPRECISION
#ifdef DOUBLEPRECISION
#define FFTW(f) fftw_ ## f
#else
#define FFTW(f) fftwf_ ## f
#endif

enum GridMode {grid_real_space, grid_fourier_space, grid_unknown_mode};

class Grid {
 public:
  Grid(const int nc);
  ~Grid();

  void fft();
  void clear();
  
  Float* fx;
  Complex* fk;
  const size_t nc;
  GridMode mode;
 private:
  const size_t ncz;
  FFTW(plan) plan;
};

#endif
