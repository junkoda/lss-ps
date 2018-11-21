#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "config.h"
#include "lognormal.h"

void generate_delta_k(Grid* const grid,
		      const unsigned long seed,
		      const bool fix_amplitude)
{
  // Convert P(k) grid to random delta(k) grid such that
  // <|delta(k)|^2> = P(k)
  //
  // input:  grid as P(k)
  // output: grid as delta(k)

  const int nc= grid->nc;
  const size_t nckz= nc/2 + 1;
  const double boxsize= grid->boxsize;
  const double vol= boxsize*boxsize*boxsize;
  fftw_complex* const fk= (fftw_complex*) grid->fx;
  
  // P(k) = 1/V <delta(k) delta^*(k)>
  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  size_t negative= 0;
  double P_min= 0.0;

  
  for(int ix=0; ix<nc; ++ix) {
   if(2*ix == nc) continue;
   
   for(int iy=0; iy<nc; ++iy) {
    if(2*iy == nc) continue;
    
    int iz0= (ix == 0 && iy == 0);
    for(int iz=iz0; iz<nc/2; ++iz) {
      size_t index= (ix*nc + iy)*nckz + iz;

      double ampl=1.0;
      if(!fix_amplitude) {
	do
	  ampl = gsl_rng_uniform(rng);
	while(ampl == 0.0);
    
	ampl= -log(ampl);
      }

      double phase= gsl_rng_uniform(rng)*2*M_PI;
      
      double delta2= vol*fk[index][0];

      double delta_k_mag= 0.0;
      if(fk[index][0] < P_min)
	P_min= fk[index][0];
	  
      if( fk[index][0] > 0.0)
	delta_k_mag= sqrt(ampl*delta2);
      else
	negative++;

      fk[index][0]= delta_k_mag*cos(phase);
      fk[index][1]= delta_k_mag*sin(phase);
	  
    }
   }
  }

  gsl_rng_free(rng);

  //fprintf(stderr, "P_min= %e\n", P_min);
  //fprintf(stderr, "negative P(k): %zu\n", negative);

  //
  // reality condition delta(-k) = delta(k)^*
  //
  for(int ix=0; ix<nc; ++ix) {
    if(2*ix == nc) continue;
    int iix= ix == 0 ? 0 :  nc - ix;
    assert(0 <= iix && iix < nc);
    for(int iy=0; iy<nc; ++iy) {
      if(2*iy == nc) continue;
      int iiy= iy == 0 ? 0 : nc - iy;

      size_t index= (ix*nc + iy)*nckz;    // index of k
      size_t iindex= (iix*nc + iiy)*nckz; // index of -k

      fk[iindex][0]= fk[index][0];
      fk[iindex][1]= -fk[index][1];
    }
  }
}

