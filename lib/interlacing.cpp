#include <cmath>
#include <complex>

#include "config.h"
#include "grid.h"

using namespace std;

void interlacing(GridComplex* const grid)
{
  const size_t nc= grid->nc;
  const int nq= nc/2 + 1;
  complex<Float>* const d= (complex<Float>*) grid->fk;

  const double tpiN= 2.0*M_PI/nc;
  const double piN= -M_PI/nc;
  complex<Float> rec= complex<Float>(cos(piN), sin(piN));

  complex<Float> c1= complex<Float>(1, 0);
  complex<Float> ci= complex<Float>(0, 1);

  complex<Float> xrec= 1;
  for(size_t ix=0; ix<nq; ++ix) {
   size_t icx= nc - ix - 1; // need to double check

   complex<Float> yrec= 1;
   for(size_t iy=0; iy<nq; ++iy) {
    size_t icy= nc - iy - 1;

    complex<Float> zrec= c1;
    for(size_t iz=0; iz<nq; ++iz) {
      size_t icz= nc - iz - 1;

      complex<Float> cma= ci*xrec*yrec*zrec;
      complex<Float> cmb= ci*xrec*yrec*conj(zrec);
      complex<Float> cmc= ci*xrec*conj(yrec)*zrec;
      complex<Float> cmd= ci*xrec*conj(yrec*zrec);
                   
      complex<Float> c000= d[(ix*nc + iy)*nc + iz]*(c1 - cma)
	                   + conj(d[(icx*nc +  icy)*nc + icz])*(c1 + cma);
      complex<Float> c001= d[(ix*nc + iy)*nc + icz]*(c1 - cmb)
	                   + conj(d[(icx*nc + icy)*nc + iz])*(c1 + cmb);
      complex<Float> c010= d[(ix*nc + icy)*nc + iz]*(c1 - cmc)
	                   + conj(d[(icx*nc + iy)*nc + icz])*(c1 + cmc);
      complex<Float> c011= d[(ix*nc + icy)*nc + icz]*(c1 - cmd)
	                   + conj(d[(icx + iy)*nc + iz])*(c1 + cmd);
                 
      d[( ix*nc + iy )*nc + iz ]= c000;
      d[( ix*nc + iy )*nc + icz]= c001;
      d[( ix*nc + icy)*nc + iz ]= c010;
      d[( ix*nc + icy)*nc + icz]= c011;
      d[(icx*nc + iy )*nc + iz ]= conj(d[(ix*nc + icy)*nc + icz]);
      d[(icx*nc + iy )*nc + icz]= conj(d[(ix*nc + icy)*nc + icz]);
      d[(icx*nc + icy)*nc + iz ]= conj(d[(ix*nc + iy )*nc + icz]);
      d[(icx*nc + icy)*nc + icz]= conj(d[(ix*nc + iy )*nc + iz ]);

      zrec *= rec;
    }
    yrec *= rec;
   }
   xrec *= rec;
  }

}
