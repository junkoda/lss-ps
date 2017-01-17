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

  for(Catalogue::const_iterator
	p= cat->begin(); p != cat->end(); ++p) {
    for(int j=0; j<3; ++j) {
      ix[j]= (int) floor((p->x[j] - x0[j])*dx_inv);
      w[j]= 1 - (p->x[j]*dx_inv - ix[j]);   // CIC weight for left point
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
	     "# density total %le %le; rel difference %e\n", total, np, err);
  assert(err < 1.0e-5);
}

