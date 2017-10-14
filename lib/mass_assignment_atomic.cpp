//
// Serial mass assignment
// Only for measuring parallelisation performance
//



#include <iostream>
#include <vector>
#include "catalogue.h"
#include "grid.h"

namespace {
//
// Mass assignment functions
//  x: position relative to the cubic box corner in units of grid spacing;
//     i.e., (0, 0, 0) and (nc, nc, nc) is the edge of the FFT grid,
//     where nc is the number of grids per dimension

struct NGP {
  void operator()(const double x[], const double w, Grid* const d) const {
    int ix[3];
    for(int i=0; i<3; ++i)
      ix[i] = (int) floor(x[i] + 0.5);

    d->add_atomic(ix[0], ix[1], ix[2], w);
  }
  static const int n_mas = 1;
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


    d->add_atomic(ix0[0], ix0[1], ix0[2], w*w0[0]*w0[1]*w0[2]);
    d->add_atomic(ix0[0], ix1[1], ix0[2], w*w0[0]*w1[1]*w0[2]);
    d->add_atomic(ix0[0], ix0[1], ix1[2], w*w0[0]*w0[1]*w1[2]);
    d->add_atomic(ix0[0], ix1[1], ix1[2], w*w0[0]*w1[1]*w1[2]);

    d->add_atomic(ix1[0], ix0[1], ix0[2], w*w1[0]*w0[1]*w0[2]);
    d->add_atomic(ix1[0], ix1[1], ix0[2], w*w1[0]*w1[1]*w0[2]);
    d->add_atomic(ix1[0], ix0[1], ix1[2], w*w1[0]*w0[1]*w1[2]);
    d->add_atomic(ix1[0], ix1[1], ix1[2], w*w1[0]*w1[1]*w1[2]);
  }

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
      for(int diy=0; diy<3; ++diy) {
	int iy= (ix0[1] + diy - 1 + d->nc) % d->nc;
	for(int diz=0; diz<3; ++diz) {
	  int iz= (ix0[2] + diz - 1 + d->nc) % d->nc;
	  d->add_atomic(ix, iy, iz, w*ww[0][dix]*ww[1][diy]*ww[2][diz]);
	}
      }
    }
  }

  static const int n_mas = 3;
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
  const double boxsize= grid->boxsize;
  const double dx_inv= nc/boxsize;
  const double dx= boxsize/nc;
  const double x0[]= {grid->x0_box[0] + grid->offset*dx,
		      grid->x0_box[1] + grid->offset*dx,
		      grid->x0_box[2] + grid->offset*dx};

  long double w_sum = 0.0;
  long double w2_sum = 0.0;
  long double nw2_sum = 0.0;

  #pragma omp parallel for
  for(size_t i=0; i<np; ++i) {
    Float rx[3];
    double w = weight == nullptr ? 1.0 : *weight;
    double nb = nbar == nullptr ? 1.0 : *nbar;
    
    const double w2 = w*w;
    w_sum += w;
    w2_sum += w2;
    nw2_sum += nb*w2;
    
    rx[0] = (xyz[0] - x0[0])*dx_inv;
    rx[1] = (xyz[1] - x0[1])*dx_inv;
    rx[2] = (xyz[2] - x0[2])*dx_inv;
    
    f(rx, w, grid);
      
    xyz    = (float_type*) ((char*) xyz    + xyz_stride);
    
    if(weight)
      weight = (float_type*) ((char*) weight + weight_stride);
    if(nbar)
      nbar   = (float_type*) ((char*) nbar   + nbar_stride);
  }
  
  grid->total_weight = w_sum;
  grid->w2_sum = w2_sum;
  grid->nw2_sum = nw2_sum;
  grid->np = np;
  grid->n_mas = f.n_mas;
}

} // end of unnamed namespace

//
// Wrapping mass_assignment_template with variable mas
//
void mass_assignment_from_array_atomic(double const * const xyz,
				const size_t xyz_stride,
				double const * const weight,
				const size_t weight_stride,
				double const * const nbar,
				const size_t nbar_stride,
				const size_t np,
				const int mas,
				bool parallelise,
				Grid* const grid)
{
  // parallelise = true if you want to parallelise the mass assignment
  // within this function call

  std::cerr << "mass_assignment_serial\n";

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



