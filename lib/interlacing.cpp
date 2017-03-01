#include <iostream>
#include <complex>
#include <vector>
#include <chrono>
#include <cmath>
#include <cassert>


#include "msg.h"
#include "config.h"
#include "grid.h"

using namespace std;

static vector<complex<Float>> exp_kh;

// Interlacing alias correction using two real grids
// Input:
//    grid: grid with offset = 0 in Fourier space 
//    grid_shifted: grid with offset = 0.5 in Fourier space
// Output:
//    grid is modified with anti-aliased Fourier modes


void interlacing_init(const int nc)
{
  if(exp_kh.size() == nc)
    return;
  
  const Float fac= M_PI/nc;
  
  for(int i=0; i<nc; ++i) {
    Float ik= i <= nc/2 ? i : i - nc;
    exp_kh.push_back(complex<Float>(cos(ik*fac), sin(ik*fac)));
  }
}

void interlacing(Grid* const grid, Grid const * const grid_shifted)
{
  // Perform interlacing
  //   Sefusatti et al, 2016, MNRAS 460 3624
  //   Accurate estimators of correlation functions in Fourier space
  //
  // Input:
  //   grid: delta[k] in Fourier space
  //   grid_shifted: shifted grid in Fourier space x -> x + 0.5*boxsize/nc
  // Output:
  //   grid: delta[k] with interlacing -- aliasing is reduced
  //
  auto ts = std::chrono::high_resolution_clock::now();
  
  assert(grid->nc == grid_shifted->nc);
  assert(grid->mode == grid_fourier_space);
  assert(grid_shifted->mode == grid_fourier_space);

  const int nc= grid->nc;
  // exp_kh = exp[(kx + ky + kz)*h] where
  // k_i = 2pi/L * ik_i
  // h = 0.5*boxsize/nc

  
  interlacing_init(nc);
  assert(exp_kh.size() == nc);

  const size_t nckz= nc/2 + 1;
  complex<Float>* const d= (complex<Float>*) grid->fk;
  complex<Float> const * const ds= (complex<Float> const*) grid_shifted->fk;

  for(int ix=0; ix<nc; ++ix) {
   complex<Float> exp_x= exp_kh[ix];
   
   for(int iy=0; iy<nc; ++iy) {
    complex<Float> exp_y= exp_kh[iy];
    complex<Float> exp_xy = exp_x*exp_y;

    for(int iz=0; iz<nckz; ++iz) {
      complex<Float> expkh = exp_xy*exp_kh[iz];

      size_t index= nckz*(nc*ix + iy) + iz;
      d[index] = 0.5*(d[index] + ds[index]*expkh);
    }
   }
  }

  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time interlacing %le\n",
	     std::chrono::duration<double>(te - ts).count());
}
