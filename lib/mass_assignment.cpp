#include <cmath>
#include <cassert>

#include "msg.h"
#include "catalogue.h"
#include "grid.h"

void mass_assignment_cic(Catalogue const * const cat,
			 const Float x0[], const Float boxsize,
			 Grid* const grid)
{
  msg_printf(msg_verbose, "Mass assigment with CIC\n");
  grid->clear();
  
  const size_t nc = grid->nc;
  const Float dx_inv= nc/boxsize;
  const size_t ncz= 2*(nc/2+1);
  Float* const d= grid->fx;
  grid->boxsize= boxsize;

  int ix[3];
  size_t ix0[3], ix1[3];
  Float w[3];
  Float rx[3];
  
  for(Catalogue::const_iterator
	p= cat->begin(); p != cat->end(); ++p) {

    for(int j=0; j<3; ++j) {
      rx[j]= (p->x[j] - x0[j])*dx_inv;
      ix[j]= (int) floor(rx[j]);
      w[j]= 1 - (rx[j] - ix[j]);            // CIC weight for left point
      ix0[j]= (ix[j] + nc) % nc;            // left grid (periodic)
      ix1[j]= (ix[j] + 1 + nc) % nc;        // right grid (periodic)
    }

    d[(ix0[0]*nc + ix0[1])*ncz + ix0[2]] += w[0]*w[1]*w[2];
    d[(ix0[0]*nc + ix1[1])*ncz + ix0[2]] += w[0]*(1-w[1])*w[2];
    d[(ix0[0]*nc + ix0[1])*ncz + ix1[2]] += w[0]*w[1]*(1-w[2]);
    d[(ix0[0]*nc + ix1[1])*ncz + ix1[2]] += w[0]*(1-w[1])*(1-w[2]);

    d[(ix1[0]*nc + ix0[1])*ncz + ix0[2]] += (1-w[0])*w[1]*w[2];
    d[(ix1[0]*nc + ix1[1])*ncz + ix0[2]] += (1-w[0])*(1-w[1])*w[2];
    d[(ix1[0]*nc + ix0[1])*ncz + ix1[2]] += (1-w[0])*w[1]*(1-w[2]);
    d[(ix1[0]*nc + ix1[1])*ncz + ix1[2]] += (1-w[0])*(1-w[1])*(1-w[2]);
  }

  // check total & normalize
  const double np= cat->size();

  double total= 0.0;
  float nbar_inv= (double)nc*nc*nc/np;

  for(int ix=0; ix<nc; ++ix) {
    for(int iy=0; iy<nc; ++iy) {
      for(int iz=0; iz<nc; ++iz) {
	size_t index= ncz*(nc*ix + iy) + iz;
	total += d[index];

  	d[index]= d[index]*nbar_inv - 1.0f; 
      }
    }
  }

  grid->shot_noise= boxsize*boxsize*boxsize/np;

  float err= fabs(total - np)/np;
  msg_printf(msg_debug,
	     "Density total %le, expected %le; rel difference %e\n",
	     total, np, err);
  assert(err < 1.0e-5);
}

//
// Interlacing
//
void mass_assignment_interlacing_cic(Catalogue const * const cat,
			 const Float x0[], const Float boxsize,
			 GridComplex* const grid)
{
  msg_printf(msg_verbose, "Mass assigment with CIC (interlacing)\n");
  grid->clear();
  
  const size_t nc = grid->nc;
  const Float dx_inv= nc/boxsize;
  Complex* const d= grid->fx;
  grid->boxsize= boxsize;

  int ix[3];
  size_t ix0[3], ix1[3];
  Float w[3];

  for(Catalogue::const_iterator
	p= cat->begin(); p != cat->end(); ++p) {

    for(int j=0; j<3; ++j) {
      ix[j]= (int) floor((p->x[j] - x0[j])*dx_inv);
      w[j]= 1 - (p->x[j]*dx_inv - ix[j]);   // CIC weight for left point
      ix0[j]= (ix[j] + nc) % nc;            // left grid (periodic)
      ix1[j]= (ix[j] + 1 + nc) % nc;        // right grid (periodic)
    }

    d[(ix0[0]*nc + ix0[1])*nc + ix0[2]][0] += w[0]*w[1]*w[2];
    d[(ix0[0]*nc + ix1[1])*nc + ix0[2]][0] += w[0]*(1-w[1])*w[2];
    d[(ix0[0]*nc + ix0[1])*nc + ix1[2]][0] += w[0]*w[1]*(1-w[2]);
    d[(ix0[0]*nc + ix1[1])*nc + ix1[2]][0] += w[0]*(1-w[1])*(1-w[2]);

    d[(ix1[0]*nc + ix0[1])*nc + ix0[2]][0] += (1-w[0])*w[1]*w[2];
    d[(ix1[0]*nc + ix1[1])*nc + ix0[2]][0] += (1-w[0])*(1-w[1])*w[2];
    d[(ix1[0]*nc + ix0[1])*nc + ix1[2]][0] += (1-w[0])*w[1]*(1-w[2]);
    d[(ix1[0]*nc + ix1[1])*nc + ix1[2]][0] += (1-w[0])*(1-w[1])*(1-w[2]);

    // Shifted grid
    for(int j=0; j<3; ++j) {
      ix[j]= (int) floor((p->x[j] - x0[j] + 0.5)*dx_inv);
      w[j]= 1 - (p->x[j]*dx_inv - ix[j]);   // CIC weight for left point
      ix0[j]= (ix[j] + nc) % nc;            // left grid (periodic)
      ix1[j]= (ix[j] + 1 + nc) % nc;        // right grid (periodic)
    }

    d[(ix0[0]*nc + ix0[1])*nc + ix0[2]][1] += w[0]*w[1]*w[2];
    d[(ix0[0]*nc + ix1[1])*nc + ix0[2]][1] += w[0]*(1-w[1])*w[2];
    d[(ix0[0]*nc + ix0[1])*nc + ix1[2]][1] += w[0]*w[1]*(1-w[2]);
    d[(ix0[0]*nc + ix1[1])*nc + ix1[2]][1] += w[0]*(1-w[1])*(1-w[2]);

    d[(ix1[0]*nc + ix0[1])*nc + ix0[2]][1] += (1-w[0])*w[1]*w[2];
    d[(ix1[0]*nc + ix1[1])*nc + ix0[2]][1] += (1-w[0])*(1-w[1])*w[2];
    d[(ix1[0]*nc + ix0[1])*nc + ix1[2]][1] += (1-w[0])*w[1]*(1-w[2]);
    d[(ix1[0]*nc + ix1[1])*nc + ix1[2]][1] += (1-w[0])*(1-w[1])*(1-w[2]);
  }

  // check total & normalize
  const double np= cat->size();

  double total_re= 0.0, total_im= 0.0;
  float nbar_inv= (double)nc*nc*nc/np;

  for(int ix=0; ix<nc; ++ix) {
    for(int iy=0; iy<nc; ++iy) {
      for(int iz=0; iz<nc; ++iz) {
	size_t index= nc*(nc*ix + iy) + iz;
	total_re += d[index][0];
	total_im += d[index][1];
	
  	d[index][0]= d[index][0]*nbar_inv - 1;
	d[index][1]= d[index][1]*nbar_inv - 1; 
      }
    }
  }

  grid->shot_noise= boxsize*boxsize*boxsize/np;

  float err= fabs(total_re - np)/np + fabs(total_im - np)/np;
  msg_printf(msg_debug,
	     "# density total %le %le; rel difference %e\n",
	     0.5*(total_re + total_im), np, err);
  assert(err < 1.0e-5);
}

