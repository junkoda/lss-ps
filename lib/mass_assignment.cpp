#include <chrono>
#include <cmath>
#include <cassert>

#include "msg.h"
#include "catalogue.h"
#include "grid.h"
#include "mass_assignment.h"

using namespace std;


//
// Mass assignment algorithm
//

/*
template<typename MAS>
void mass_assignment_template(const vector<Particle>& cat,
			      const double x0[],
			      const double boxsize,
			      MAS f,
			      bool parallel,
			      Grid* grid)
{
  // Add density to grid
  const double dx_inv = grid->nc/boxsize;
  double w_sum = 0.0;
  double w2_sum = 0.0;
  double nw2_sum = 0.0;

#pragma omp parallel for if(parallel) reduction( + : w_sum, w2_sum, nw2_sum )
  for(size_t i=0; i<cat.size(); ++i) {
    Float rx[3];
    double w = cat[i].w;
    const double nbar = cat[i].nbar;
    const double w2 = w*w;
    w_sum += w;
    w2_sum += w2;
    nw2_sum += nbar*w2;
    
    rx[0] = (cat[i].x[0] - x0[0])*dx_inv;
    rx[1] = (cat[i].x[1] - x0[1])*dx_inv;
    rx[2] = (cat[i].x[2] - x0[2])*dx_inv;

    f(rx, w, grid);
  }
 
  #pragma omp critical (__MASS_ASSIGNMENT__)
  {
    grid->total_weight += w_sum;
    grid->w2_sum += w2_sum;
    grid->nw2_sum += nw2_sum;
    grid->np += cat.size();
    
    grid->n_mas = f.n_mas;
    grid->boxsize= boxsize;
  }
}
*/

//
// Wrapper for for vector<Particle>
//
void mass_assignment_from_particles(const vector<Particle>& cat,
				    const int mas, const bool parallelise,
				    Grid* const grid)
{
  // parallelise = true if you want to parallelise the mass assignment
  // within this function call

  Float const * const xyz= cat[0].x;
  Float const * const weight= &cat[0].w;
  Float const * const nbar= &cat[0].nbar;
  const size_t stride= sizeof(Particle);
  
  switch(mas) {
  case 1:
    //mass_assignment_template(cat, x0, boxsize, NGP(), parallelise, grid);
    mass_assignment_template(xyz, stride,
			     weight, stride,
			     nbar, stride,
			     cat.size(),
			     NGP(), parallelise, grid);
    break;

  case 2:
    mass_assignment_template(xyz, stride,
			     weight, stride,
			     nbar, stride,
			     cat.size(),
			     CIC(), parallelise, grid);
    break;

  case 3:
    mass_assignment_template(xyz, stride,
			     weight, stride,
			     nbar, stride,
			     cat.size(),
			     TSC(), parallelise, grid);
    break;

  default:
    msg_printf(msg_error, "Error: unknown mass assignment scheme %d\n",
	       mas);
  }
}


