#ifndef MOMENT_H
#define MOMENT_H 1

#include <complex>
#include <cassert>
#include <cmath>
#include "grid.h"

//
// Moment function objects
//

// (x_i/r)
struct Moment1 {
  Moment1(const int i_) : i(i_) {}

  double operator()(const double x[]) const {
    double r= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    return r > 0.0 ? x[i]/r : 0.0;
  }
  const int i;
};

// (x_i/r) (x_j/r)
struct Moment2 {
  Moment2(const int i_, const int j_) : i(i_), j(j_) {}

  double operator()(const double x[]) const {
    double r2= x[0]*x[0] + x[1]*x[1] + x[2]*x[2];

    return r2 > 0.0 ? x[i]*x[j]/r2 : 0.0;
  }
  const int i, j;
};

// (x_i/r) (x_j/r) (x_k/r)
struct Moment3 {
  Moment3(const int i_, const int j_, const int k_) : i(i_), j(j_), k(k_) {}

  double operator()(const double x[]) const {
    double r2= x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r= sqrt(r2);

    return r2 > 0.0 ? x[i]*x[j]*x[k]/(r*r2) : 0.0;
  }
  const int i, j, k;
};


// (x_i/r) (x_j/r) (x_k/r) (x_l/r)
struct Moment4 {
  Moment4(const int i_, const int j_, const int k_, const int l_) :
    i(i_), j(j_), k(k_), l(l_) {}

  double operator()(const double x[]) const {
    double r2= x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    return r2 > 0.0 ? (x[i]*x[j]/r2)*(x[k]*x[l]/r2) : 0.0;
  }
  const int i, j, k, l;
};


//
// Moment function objects for window function Qnl
//

// r^-n
struct WindowMoment0 {
  WindowMoment0(const int n_) : n(n_) {}

  double operator()(const double x[]) const {
    double r= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    return r > 0.0 ? 1.0/pow(r, n + 1) : 0.0;
  }
  const int n;
};

// (x_i/r) r^-n
struct WindowMoment1 {
  WindowMoment1(const int i_, const int n_) : i(i_), n(n_) {}

  double operator()(const double x[]) const {
    double r= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    return r > 0.0 ? x[i]/pow(r, n + 1) : 0.0;
  }
  const int i, n;
};

// (x_i/r) (x_j/r) r^-n
struct WindowMoment2 {
  WindowMoment2(const int i_, const int j_, const int n_) :
    i(i_), j(j_), n(n_) {}

  double operator()(const double x[]) const {
    double r= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    return r > 0.0 ? x[i]*x[j]/pow(r, n + 2) : 0.0;
  }
  const int i, j, n;
};

// (x_i/r) (x_j/r) (x_k/r) r^-n
struct WindowMoment3 {
  WindowMoment3(const int i_, const int j_, const int k_, const int n_) :
    i(i_), j(j_), k(k_), n(n_) {}

  double operator()(const double x[]) const {
    double r= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    return r > 0.0 ? x[i]*x[j]*x[k]/pow(r, n + 3) : 0.0;
  }
  const int i, j, k, n;
};


// (x_i/r) (x_j/r) (x_k/r) (x_l/r) r^-n
struct WindowMoment4 {
  WindowMoment4(const int i_, const int j_, const int k_, const int l_,
		const int n_): i(i_), j(j_), k(k_), l(l_), n(n_) {}

  double operator()(const double x[]) const {
    double r= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    return r > 0.0 ? x[i]*x[j]*x[k]*x[l]/pow(r, n + 4) : 0.0;
  }
  const int i, j, k, l, n;
};


//
// Compute moment
// Q_ij...(x) = x_i/r x_j/r ... delta(x)
// for Yamamoto FFT estimator
//

template<typename Moment>
void moment_compute_x(Grid const * const grid,
		      Float const * const x0,
		      Moment moment,
		      Grid* const grid_moment)
{
  // Compute moment
  //   grid_moment(x) = moment(x) grid(x),
  // Args:
  //   grid:        delta(x) grid
  //   moment:      function object Moment2 or Moment4
  //   grid_moment: output grid of moment
  assert(grid->mode == GridMode::real_space);
  
  const size_t nc= grid->nc;
  const size_t ncz= grid->ncz;
  const double dx= grid->boxsize / nc;

  double* const m= grid_moment->fx;
  double const * const d= grid->fx;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(size_t ix=0; ix<nc; ++ix) {
    double x[3];
    
    // position of the grid point
    x[0]= x0[0] + dx*ix;
    for(size_t iy=0; iy<nc; ++iy) {
      x[1]= x0[1] + dx*iy;
      for(size_t iz=0; iz<nc; ++iz) {
	x[2]= x0[2] + dx*iz;
	size_t index= (ix*nc + iy)*ncz + iz;

	m[index] = moment(x)*d[index];
      }
    }
  }

  grid_moment->mode = GridMode::real_space;
}


template<typename Moment>
void moment_compute_k(Grid const * const grid,
		      const Moment moment,
		      const double coef,
		      Grid* const grid_out)
{
  // Add moment of grid to grid_out in Fourier space
  // using function object moment
  // grid_moment(k) += coef moment(k) grid(k)
  // moment(k) = (k_i/k) (k_j/k) ...
  const int nc= static_cast<int>(grid->nc);

  assert(grid_out->nc == grid->nc);
  assert(grid->mode == GridMode::fourier_space);

  const size_t nckz= nc/2+1;
  const int iknq= nc/2;

  std::complex<double> const * const d_k= grid->fk;
  std::complex<double>* const m_k= grid_out->fk;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ix=0; ix<nc; ++ix) {
    double k[3];
    k[0] = ix <= iknq ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      k[1] = iy <= iknq ? iy : iy - nc;
      int kz0 = !(k[0] > 0.0 || (k[0] <= 0.0 && k[1] > 0.0));

      for(int iz=kz0; iz<iknq; ++iz) {
	k[2]= iz;
	size_t index= nckz*(nc*ix + iy) + iz;
	const double m= coef*moment(k);
	
	m_k[index] += m*d_k[index];	
      }
    }
  }
  grid_out->mode= GridMode::fourier_space;
}



#endif
