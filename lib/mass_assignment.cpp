#include <chrono>
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
  auto ts = std::chrono::high_resolution_clock::now();
  
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

  double w_sum = 0.0;
  double w2_sum = 0.0;
  double nw2_sum = 0.0;

  for(Catalogue::const_iterator
	p= cat->begin(); p != cat->end(); ++p) {

    double ww = p->w;
    const double nbar = p->nbar;
    
    const double w2 = ww*ww;
    w_sum += ww;
    w2_sum += w2;    
    nw2_sum += nbar*w2;

    for(int j=0; j<3; ++j) {
      rx[j]= (p->x[j] - x0[j])*dx_inv;
      ix[j]= (int) floor(rx[j]);
      w[j]= 1 - (rx[j] - ix[j]);            // CIC weight for left point
      ix0[j]= (ix[j] + nc) % nc;            // left grid (periodic)
      ix1[j]= (ix[j] + 1 + nc) % nc;        // right grid (periodic)
    }

    d[(ix0[0]*nc + ix0[1])*ncz + ix0[2]] += ww*w[0]*w[1]*w[2];
    d[(ix0[0]*nc + ix1[1])*ncz + ix0[2]] += ww*w[0]*(1-w[1])*w[2];
    d[(ix0[0]*nc + ix0[1])*ncz + ix1[2]] += ww*w[0]*w[1]*(1-w[2]);
    d[(ix0[0]*nc + ix1[1])*ncz + ix1[2]] += ww*w[0]*(1-w[1])*(1-w[2]);

    d[(ix1[0]*nc + ix0[1])*ncz + ix0[2]] += ww*(1-w[0])*w[1]*w[2];
    d[(ix1[0]*nc + ix1[1])*ncz + ix0[2]] += ww*(1-w[0])*(1-w[1])*w[2];
    d[(ix1[0]*nc + ix0[1])*ncz + ix1[2]] += ww*(1-w[0])*w[1]*(1-w[2]);
    d[(ix1[0]*nc + ix1[1])*ncz + ix1[2]] += ww*(1-w[0])*(1-w[1])*(1-w[2]);
  }

  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time mass_assignment CIC %le\n",
	     std::chrono::duration<double>(te - ts).count());
  ts = te;

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

  te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time get fluctuation %le\n",
	     std::chrono::duration<double>(te - ts).count());

  
  float err= fabs(total - np)/np;
  msg_printf(msg_debug,
	     "Density total %le, expected %le; rel difference %e\n",
	     total, np, err);
  assert(err < 1.0e-5);
}

/*

template<typename MAS>
void mass_assignment(Catalogue const * const cat,
		     MAS f, const double x0[], const double boxsize,
		     const bool useFKP, const double Pest,
		     Grid* const grid)
{
  // Assign a moment of Galaxies in vector v to grid,
  // using mass assignment function f
  
  // get start time
  auto ts = std::chrono::high_resolution_clock::now();

  msg_printf(msg_verbose,
	     "Assigning density field on grid with Teplate n_mas=%d",
	     f.n_mas);

  const double dx_inv = grid->nc/boxsize;

  double* const d = grid->fx;

  double w_sum = 0.0;
  double w2_sum = 0.0;
  double nw2_sum = 0.0;

  for(std::vector<Particle>::const_iterator p= cat->begin();
      p != cat->end(); ++p) {
    double w = p->w;
    const double nbar = 1.0; // TODO nbar
    
    if(useFKP)
       w /=  (1.0 + nbar*Pest);

    const double w2 = w*w;
    w_sum += w;
    w2_sum += w*w;    
    nw2_sum += nbar*w2;
    
    double x[] = {(p->x[0] - x0[0])*dx_inv,
		  (p->x[1] - x0[1])*dx_inv,
		  (p->x[2] - x0[2])*dx_inv};
    f(x);
  }

  grid->total_weight = w_sum;
  grid->raw_noise = w2_sum;
  grid->normalisation = nw2_sum;
  grid->n_mas = f.n_mas;
  
  // time duration
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time mass_assignment template %le\n",
	     std::chrono::duration<double>(te - ts).count());
}
*/

//}
