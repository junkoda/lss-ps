#include <iostream>
#include <complex>
#include <cmath>
#include <cassert>


#include "config.h"
#include "grid.h"

using namespace std;


// Interlacing alias correction using two real grids
// Input:
//    grid: grid with offset = 0 in Fourier space 
//    grid_shifted: grid with offset = 0.5 in Fourier space
// Output:
//    grid is modified with anti-aliased Fourier modes

void interlacing2(Grid* const grid, Grid const * const grid_shifted)
{
  assert(grid->nc == grid_shifted->nc);
  assert(grid->mode == grid_fourier_space);
  assert(grid_shifted->mode == grid_fourier_space);

  const int nc= grid->nc;
  const Float fac= 2.0*M_PI/grid->boxsize;
  const Float h= 0.5*grid->boxsize/nc;

  // Todo: boxsize here can be canceled
  const size_t nckz= nc/2 + 1;
  complex<Float>* const d= (complex<Float>*) grid->fk;
  complex<Float> const * const ds= (complex<Float> const*) grid_shifted->fk;

  for(int ix=0; ix<nc; ++ix) {
   Float kx= ix <= nc/2 ? fac*ix : fac*(ix - nc);
   
   for(int iy=0; iy<nc; ++iy) {
    Float ky= iy <= nc/2 ? fac*iy : fac*(iy - nc);
    int iz0 = !(kx > 0.0 || (kx == 0.0 && ky > 0.0));
    
    for(int iz=iz0; iz<nckz; ++iz) {
      Float kz= fac*iz;
      Float kh= (kx + ky + kz)*h;
      complex<Float> expkh = complex<Float>(cos(kh), sin(kh));

      size_t index= nckz*(nc*ix + iy) + iz;
      d[index] = 0.5*(d[index] + ds[index]*expkh);
    }
   }
  }
}
