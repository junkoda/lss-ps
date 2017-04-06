#include <chrono>
#include <cmath>
#include <cassert>

#include "msg.h"
#include "catalogue.h"
#include "grid.h"

using namespace std;

//
// Mass assignment functions
//  x: position relative to the cubic box corner in units of grid spacing;
//     i.e., (0, 0, 0) and (nc, nc, nc) is the edge of the FFT grid,
//     where nc is the number of grids per dimension

struct NGP {
  void operator()(const double x[], const double w, Grid* const d) const {
    int ix[3];
    for(int k=0; k<3; ++k) {
      ix[k] = (int)(x[k] + 0.5);
    }

    d->add(ix[0], ix[1], ix[2], w);
  }
  static const int n_mas = 1;
};
    

struct CIC {
  void operator()(const double x[], const double w, Grid* const d) const {
    int ix[3], ix0[3], ix1[3];
    double w0[3], w1[3];
    for(int k=0; k<3; ++k) {
      ix[k] = (int) x[k]; // assume x[k] >= 0.0
      ix0[k]= (ix[k] + d->nc) % d->nc;    // left grid point (periodic)
      ix1[k]= (ix[k] + 1 + d->nc) % d->nc;// right grid point (periodic)

      w1[k] = x[k] - ix[k];              // CIC weight to right grid point
      w0[k] = 1.0 - w1[k];                 //               left grid point
    }
    
    d->add(ix0[0], ix0[1], ix0[2], w*w0[0]*w0[1]*w0[2]);
    d->add(ix0[0], ix1[1], ix0[2], w*w0[0]*w1[1]*w0[2]);
    d->add(ix0[0], ix0[1], ix1[2], w*w0[0]*w0[1]*w1[2]);
    d->add(ix0[0], ix1[1], ix1[2], w*w0[0]*w1[1]*w1[2]);

    d->add(ix1[0], ix0[1], ix0[2], w*w1[0]*w0[1]*w0[2]);
    d->add(ix1[0], ix1[1], ix0[2], w*w1[0]*w1[1]*w0[2]);
    d->add(ix1[0], ix0[1], ix1[2], w*w1[0]*w0[1]*w1[2]);
    d->add(ix1[0], ix1[1], ix1[2], w*w1[0]*w1[1]*w1[2]);
  }

  static const int n_mas = 2;
};

struct TSC {
  void operator()(const double x[], const double w, Grid* const d) const {
    int ix0[3];

    double ww[3][3];
    for(int k=0; k<3; ++k) {
      ix0[k] = (int) (x[k] + 0.5);
      double dx1 = x[k] - ix0[k];
      double dx2 = 0.5 - dx1;
      
      ww[k][0] = 0.5*dx2*dx2;
      ww[k][1] = 0.75 - dx1*dx1;
      ww[k][2] = 1.0 - ww[k][0] - ww[k][1];
    }

    for(int dix=0; dix<3; ++dix) {
      int ix= (ix0[0] + dix - 1 + d->nc) % d->nc;
      for(int diy=0; diy<3; ++diy) {
	int iy= (ix0[1] + diy - 1 + d->nc) % d->nc;
	for(int diz=0; diz<3; ++diz) {
	  int iz= (ix0[2] + diz - 1 + d->nc) % d->nc;
	  d->add(ix, iy, iz, w*ww[0][dix]*ww[1][diy]*ww[2][diz]);
	}
      }
    }
  }

  static const int n_mas = 3;
};

//
// Mass assignment algorithm
//

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


//
// Wrapper for various mass assignment scheme
//
void mass_assignment(const vector<Particle>& cat,
		     const Float x0[], const Float boxsize,
		     const int mas, const bool parallelise,
		     Grid* const grid)
{
  // parallelise = true if you want to parallelise the mass assignment
  // within this function call
  switch(mas) {
  case 1:
    mass_assignment_template(cat, x0, boxsize, NGP(), parallelise, grid);
    break;

  case 2:
    mass_assignment_template(cat, x0, boxsize, CIC(), parallelise, grid);
    break;

  case 3:
    mass_assignment_template(cat, x0, boxsize, TSC(), parallelise, grid);
    break;

  default:
    msg_printf(msg_error, "Error: unknown mass assignment scheme %d\n",
	       mas);
  }
}

