#ifndef GRID_H
#define GRID_H 1

#include <complex>
#include <fftw3.h>
#include "config.h"

// FFTW() macro adds fftw_ or fftwf_ prefix depending on DOUBLEPRECISION
#ifdef DOUBLEPRECISION
#define FFTW(f) fftw_ ## f
#else
#define FFTW(f) fftwf_ ## f
#endif

enum class GridMode {unknown_mode, real_space, fourier_space};

//
// Class Grid
//
class Grid {
 public:
  Grid(const int nc);
  ~Grid();

  double& operator[](const size_t i) {
    return fx[i];
  }

  double operator[](const size_t i) const {
    return fx[i];
  }

  double& operator()(const size_t ix, const size_t iy, const size_t iz) {
    return fx[(ix*nc + iy)*ncz + iz];
  }

  double operator()(const size_t ix, const size_t iy, const size_t iz) const {
    return fx[(ix*nc + iy)*ncz + iz];
  }

  void add(const size_t ix, const size_t iy, const size_t iz,
	   const Float val) {
    size_t index= (ix*nc + iy)*ncz + iz;
    #pragma omp atomic
    fx[index] += val;
  }
  
  void fft();
  void clear();
  void write(const char filename[]);
  
  Float* fx;
  std::complex<Float>* fk;
  const size_t nc, ncz;
  GridMode mode;
  Float boxsize;
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

//
// Moment function objects
//
/*
struct Moment2 {
  Moment2(const int i_, const int j_) : i(i_), j(j_) {}
  Float operator()(const Float x[]) const {
    Float r2= x*x + y*y + z*z;
    return x[i]*x[j]/r2;
  }
  const int i, j;
};

struct Moment4 {
Moment3(const int i_, const int j_, const int k_, const int l_) :
  i(i_), j(j_), k(k_), l(l_) {}
  Float operator()(const Float x[]) const {
    Float r2= x*x + y*y + z*z;
    return (x[i]*x[j]/r2)*(x[k]*x[l]/r2);
  }
  const int i,j,k,l;
};
*/


//
// Compute grid of moments for Yamamoto estimator
//
/*
template<Moment moment>
void compute_moment(const Grid& grid_delta, Grid& grid_moment)
{
  const size_t nc= grid_delta.nc;
  const size_t ncz= grid_delta.ncz;
  const double dx= grid_delta.boxsize / nc;
  double const * const x0= grid_delta.x0;

  double* const m= grid_moment.fx;
  double const * const d= grid_delta.fx;

  double x[3];
  
  for(size_t ix=0; ix<nc; ++ix) {
   double x[0]= x0[0] + dx*ix;
   for(size_t iy=0; iy<nc; ++iy) {
    double x[1]= x0[1] + dx*iy;
    for(size_t iz=0; iz<nc; ++iz) {
     double x[2]= x0[2] + dx*iz;
     size_t index= (ix*nc + iy)*ncz + iz;
     m[index] = moment(x)*d[index];
    }
   }
  }
}
*/

#endif
