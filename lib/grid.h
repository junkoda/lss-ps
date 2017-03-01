#ifndef GRID_H
#define GRID_H 1

#include <fftw3.h>
#include "config.h"
#include "grid.h"

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

  double& operator[](const size_t i) {
    return fx[i];
  }

  double& operator[](const size_t i) const {
    return fx[i];
  }

  double& operator()(const size_t ix, const size_t iy, const size_t iz) {
    return fx[(ix*nc + iy)*ncz + iz];
  }

  double& operator()(const size_t ix, const size_t iy, const size_t iz) const {
    return fx[(ix*nc + iy)*ncz + iz];
  }

  void fft();
  void clear();
  void write(const char filename[]);
  
  Float* fx;
  Complex* fk;
  const size_t nc, ncz;
  GridMode mode;
  Float boxsize;
  Float x0[3];
  Float shot_noise;
  Float pk_normalisation;

  double total_weight, w2_sum, nw2_sum;
  size_t np;
  int n_mas;

 private:
  FFTW(plan) plan;  
};

void grid_print_time();

void grid_compute_fluctuation(Grid& grid1, Grid& grid_rand);
void grid_compute_fluctuation_homogeneous(Grid& grid1);


#endif
