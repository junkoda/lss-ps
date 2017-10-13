#ifndef MASS_ASSIGNMENT_H
#define MASS_ASSIGNMENT_H 1

#include <iostream>
#include <vector>
#include "catalogue.h"
#include "grid.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//
// Mass assignment functions
//  x: position relative to the cubic box corner in units of grid spacing;
//     i.e., (0, 0, 0) and (nc, nc, nc) is the edge of the FFT grid,
//     where nc is the number of grids per dimension

struct NGP {
  void operator()(const double x[], const double w, Grid* const d) const {
    int ix[3];
    for(int k=0; k<3; ++k) {
      ix[k] = (int) floor(x[k] + 0.5);
    }

    if(ix_left <= ix[0] && ix[0] < ix_right)
      d->add(ix[0], ix[1], ix[2], w);
  }
  static const int n_mas = 1;
  int ix_left, ix_right;
};
    

struct CIC {
  void operator()(const double x[], const double w, Grid* const d) const {
    int ix[3], ix0[3], ix1[3];
    double w0[3], w1[3];

    for(int k=0; k<3; ++k) {
      ix[k] = (int) floor(x[k]);
      ix0[k]= (ix[k] + d->nc) % d->nc;    // left grid point (periodic)
      ix1[k]= (ix[k] + 1 + d->nc) % d->nc;// right grid point (periodic)

      w1[k] = x[k] - ix[k];              // CIC weight to right grid point
      w0[k] = 1.0 - w1[k];                 //               left grid point
    }

    if(ix_left <= ix0[0] && ix0[0] < ix_right) {
      d->add(ix0[0], ix0[1], ix0[2], w*w0[0]*w0[1]*w0[2]);
      d->add(ix0[0], ix1[1], ix0[2], w*w0[0]*w1[1]*w0[2]);
      d->add(ix0[0], ix0[1], ix1[2], w*w0[0]*w0[1]*w1[2]);
      d->add(ix0[0], ix1[1], ix1[2], w*w0[0]*w1[1]*w1[2]);
    }

    if(ix_left <= ix1[0] && ix1[0] < ix_right) {
      d->add(ix1[0], ix0[1], ix0[2], w*w1[0]*w0[1]*w0[2]);
      d->add(ix1[0], ix1[1], ix0[2], w*w1[0]*w1[1]*w0[2]);
      d->add(ix1[0], ix0[1], ix1[2], w*w1[0]*w0[1]*w1[2]);
      d->add(ix1[0], ix1[1], ix1[2], w*w1[0]*w1[1]*w1[2]);
    }
  }

  int ix_left, ix_right;
  static const int n_mas = 2;
};

struct TSC {
  void operator()(const double x[], const double w, Grid* const d) const {
    int ix0[3];

    double ww[3][3];
    for(int k=0; k<3; ++k) {
      ix0[k] = (int) floor(x[k] + 0.5);
      double dx1 = x[k] - ix0[k];
      double dx2 = 0.5 - dx1;
      
      ww[k][0] = 0.5*dx2*dx2;
      ww[k][1] = 0.75 - dx1*dx1;
      ww[k][2] = 1.0 - ww[k][0] - ww[k][1];
    }

    for(int dix=0; dix<3; ++dix) {
      int ix= (ix0[0] + dix - 1 + d->nc) % d->nc;
      if(ix_left <= ix && ix < ix_right) {
	for(int diy=0; diy<3; ++diy) {
	  int iy= (ix0[1] + diy - 1 + d->nc) % d->nc;
	  for(int diz=0; diz<3; ++diz) {
	    int iz= (ix0[2] + diz - 1 + d->nc) % d->nc;
	    d->add(ix, iy, iz, w*ww[0][dix]*ww[1][diy]*ww[2][diz]);
	  }
	}
      }
    }
  }

  static const int n_mas = 3;
  int ix_left, ix_right;
};

//
// The mass assignment algorithm
//
template<typename MAS, typename float_type>
void mass_assignment_template(float_type const * xyz,
			      const size_t xyz_stride,
			      float_type const * weight,
			      const size_t weight_stride,
			      float_type const * nbar,
			      const size_t nbar_stride,
			      const size_t np,
			      MAS f,
			      bool parallelise,
			      Grid* const grid)
{
  // Add density to grid (grid will not be reset)
  //
  // Args:
  //     xyz:         pointer to xyz coordinate
  //     xyz_stride:  number of bytes between xyz data
  //     weight:      pointer to weight (can be null)
  //     nbar:        pointer to nbar (can be null)
  //     np:          number of particles
  //     parallelise: OpenMP paralesllise the mass assignment loop
  //     grid:        the density grid
  //
  // Output:
  //     update
  //     grid->total_weight;
  //     grid->w2_sum;
  //     grid->nw2_sum;
  //     grid->np;
  //     grid->n_mas = f.n_mas;

  assert(xyz);
  assert(grid->boxsize > 0);
  assert(grid->nc > 0);

  const int nc= grid->nc;
  const double dnc = nc;
  const double boxsize= grid->boxsize;
  const double dx_inv= nc/boxsize;
  const double dx= boxsize/nc;
  const double x0[]= {grid->x0_box[0] + grid->offset*dx,
		      grid->x0_box[1] + grid->offset*dx,
		      grid->x0_box[2] + grid->offset*dx};

  #pragma omp parallel firstprivate(xyz, weight, nbar)
  {
    int ithread= omp_get_thread_num();
    int nthread= omp_get_num_threads();

    int ix_left= ithread*nc/nthread;
    int ix_right= (ithread + 1)*nc/nthread;

    double x_left= (ix_left - 0.5*f.n_mas);
    double x_right= (ix_right + 0.5*f.n_mas);

    MAS f_local;
    f_local.ix_left = ix_left;
    f_local.ix_right= ix_right;

    #pragma omp critical
    {
      std::cerr << ithread << " / " << nthread << std::endl;
    }

    long double w_sum = 0.0;
    long double w2_sum = 0.0;
    long double nw2_sum = 0.0;

    for(size_t i=0; i<np; ++i) {
      Float rx[3];
      double w = weight == nullptr ? 1.0 : *weight;
      double nb = nbar == nullptr ? 1.0 : *nbar;
      
      const double w2 = w*w;
      w_sum += w;
      w2_sum += w2;
      nw2_sum += nb*w2;
    
      rx[0] = (xyz[0] - x0[0])*dx_inv;

      // this should not be correct due to periodic wrapup.
      // DEBUG
      //  unit of x_left is wrong!!!
      if((x_left <= rx[0] && rx[0] <= x_right) ||
	 (rx[0] >= x_left + dnc) ||
	 (rx[0] <= x_right + dnc)) {
	rx[1] = (xyz[1] - x0[1])*dx_inv;
	rx[2] = (xyz[2] - x0[2])*dx_inv;

	f_local(rx, w, grid);
      }
      
      xyz    = (float_type*) ((char*) xyz    + xyz_stride);

      if(weight)
	weight = (float_type*) ((char*) weight + weight_stride);
      if(nbar)
	nbar   = (float_type*) ((char*) nbar   + nbar_stride);
    }

    if(ithread == 0) {
      grid->total_weight = w_sum;
      grid->w2_sum = w2_sum;
      grid->nw2_sum = nw2_sum;
      grid->np = np;
      grid->n_mas = f.n_mas;
    }
  }
}


/*
void mass_assignment(const std::vector<Particle>& cat,
		     const Float x0[], const Float boxsize,
		     const int mas,
		     const bool parallel,
		     Grid* const grid);
*/
//
// Wrapper of mass_assignment_template for vector<Particle>
//
void mass_assignment_from_particles(const std::vector<Particle>& cat,
				    const int mas, const bool parallelise,
				    Grid* const grid);

//
// Wrapping mass_assignment_template with variable mas
//
template<typename float_type>
void mass_assignment_from_array(float_type const * const xyz,
				const size_t xyz_stride,
				float_type const * const weight,
				const size_t weight_stride,
				float_type const * const nbar,
				const size_t nbar_stride,
				const size_t np,
				const int mas,
				bool parallelise,
				Grid* const grid)
{
  // parallelise = true if you want to parallelise the mass assignment
  // within this function call

  //std::cerr << "mass_assignment_from_array\n";
  //std::cerr << "parallelise " << parallelise << std::endl;
  //std::cerr << "omp_get_max_threads " << omp_get_max_threads() << std::endl;

  switch(mas) {
  case 1:
    mass_assignment_template(xyz,    xyz_stride,
			     weight, weight_stride,
			     nbar,   nbar_stride,
			     np,
			     NGP(), parallelise, grid);
    break;

  case 2:
    mass_assignment_template(xyz,    xyz_stride,
			     weight, weight_stride,
			     nbar,   nbar_stride,
			     np,
			     CIC(), parallelise, grid);
    break;

  case 3:
    mass_assignment_template(xyz,    xyz_stride,
			     weight, weight_stride,
			     nbar,   nbar_stride,
			     np,
			     TSC(), parallelise, grid);
    break;

  default:
    msg_printf(msg_error, "Error: unknown mass assignment scheme %d\n",
	       mas);
    throw TypeError();
  }
}


#endif
