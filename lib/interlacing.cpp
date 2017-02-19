#include <iostream>
#include <complex>
#include <cmath>
#include <cassert>


#include "config.h"
#include "grid.h"

using namespace std;

// Intrlacing alias correction using one complex grid
void interlacing(GridComplex* const grid)
{
  const size_t nc= grid->nc;
  const size_t nq= nc/2 + 1;
  complex<Float>* const d= (complex<Float>*) grid->fk;

  const double piN= M_PI/nc;

  complex<Float> rec= complex<Float>(cos(piN), sin(piN));

  complex<Float> c1= complex<Float>(1, 0);
  complex<Float> ci= complex<Float>(0, 1);

  complex<Float> xrec= 1;
  for(size_t ix=0; ix<nq; ++ix) {
  size_t icx= (nc - ix) % nc;

   complex<Float> yrec= 1;
   for(size_t iy=0; iy<nq; ++iy) {
   size_t icy= (nc - iy) % nc;

    complex<Float> zrec= 1;
    for(size_t iz=0; iz<nq; ++iz) {
      size_t icz= (nc - iz) % nc;

      complex<Float> cma= ci*xrec*yrec*zrec;
      complex<Float> cmb= ci*xrec*yrec*conj(zrec);
      complex<Float> cmc= ci*xrec*conj(yrec)*zrec;
      complex<Float> cmd= ci*xrec*conj(yrec*zrec);
                   
      complex<Float> c000=        d[(ix *nc + iy )*nc + iz ] *(c1 - cma)
	                   + conj(d[(icx*nc + icy)*nc + icz])*(c1 + cma);
      complex<Float> c001=        d[(ix *nc + iy )*nc + icz] *(c1 - cmb)
	                   + conj(d[(icx*nc + icy)*nc + iz ])*(c1 + cmb);
      complex<Float> c010=        d[(ix *nc + icy)*nc + iz ] *(c1 - cmc)
	                   + conj(d[(icx*nc + iy )*nc + icz])*(c1 + cmc);
      complex<Float> c011=        d[(ix *nc + icy)*nc + icz] *(c1 - cmd)
	                   + conj(d[(icx*nc + iy) *nc + iz ])*(c1 + cmd);
                 
      d[( ix*nc + iy )*nc + iz ]= c000;
      d[( ix*nc + iy )*nc + icz]= c001;
      d[( ix*nc + icy)*nc + iz ]= c010;
      d[( ix*nc + icy)*nc + icz]= c011;
      d[(icx*nc + iy )*nc + iz ]= conj(d[(ix*nc + icy)*nc + icz]);
      d[(icx*nc + iy )*nc + icz]= conj(d[(ix*nc + icy)*nc + iz]);
      d[(icx*nc + icy)*nc + iz ]= conj(d[(ix*nc + iy )*nc + icz]);
      d[(icx*nc + icy)*nc + icz]= conj(d[(ix*nc + iy )*nc + iz ]);

      zrec *= rec;
    }
    yrec *= rec;
   }
   xrec *= rec;
  }

}

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
