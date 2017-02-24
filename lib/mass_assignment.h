#ifndef MASS_ASSIGNMENT_H
#define MASS_ASSIGNMENT_H 1

#include "catalogue.h"
#include "grid.h"


//
// Mass assignment functions
//  x: position relative to the cubic box corner in units of grid spacing;
//     i.e., (0, 0, 0) and (nc, nc, nc) is the edge of the FFT grid,
//     where nc is the number of grids per dimension

struct NGP {
  inline void operator()(const double x[], const double w, Grid& d) {
#ifdef USEEIGEN
    abort();
#else
    int ix[3];
    for(int k=0; k<3; ++k) {
      ix[k] = (int) round(x[k]);
    }
#endif

    d(ix[0], ix[1], ix[2]) += w;
  }
  static const int n_mas = 1;
};
    

struct CIC {
  inline void operator()(const double x[], const double w, Grid& d) {
    int ix[3], ix0[3], ix1[3];
    double w0[3], w1[3];
    for(int k=0; k<3; ++k) {
      ix[k] = (int) floor(x[k]);
      ix0[k]= (ix[k] + d.nc) % d.nc;    // left grid point (periodic)
      ix1[k]= (ix[k] + 1 + d.nc) % d.nc;// right grid point (periodic)

      w1[k] = x[k] - ix[k];              // CIC weight to right grid point
      w0[k] = 1 - w1[k];                 //               left grid point
    }
    
    d(ix0[0], ix0[1], ix0[2]) += w*w0[0]*w0[1]*w0[2];
    d(ix0[0], ix1[1], ix0[2]) += w*w0[0]*w1[1]*w0[2];
    d(ix0[0], ix0[1], ix1[2]) += w*w0[0]*w0[1]*w1[2];
    d(ix0[0], ix1[1], ix1[2]) += w*w0[0]*w1[1]*w1[2];

    d(ix1[0], ix0[1], ix0[2]) += w*w1[0]*w0[1]*w0[2];
    d(ix1[0], ix1[1], ix0[2]) += w*w1[0]*w1[1]*w0[2];
    d(ix1[0], ix0[1], ix1[2]) += w*w1[0]*w0[1]*w1[2];
    d(ix1[0], ix1[1], ix1[2]) += w*w1[0]*w1[1]*w1[2];
  }

  static const int n_mas = 2;
};

void mass_assignment_cic(Catalogue const * const cat,
			 const Float x0[], const Float boxsize,
			 const bool useFKP, const double Pest,
			 Grid* const grid);


//
// Mass assignment algorithm
//

template<typename MAS>
void mass_assignment(Catalogue const * const cat,
		     MAS f, const double x0[], const double boxsize,
		     const bool useFKP, const double Pest,
		     Grid& grid)
{
  // Assign a moment of Galaxies in vector v to grid,
  // using mass assignment function f
  
  // get start time
  auto ts = std::chrono::high_resolution_clock::now();

  msg_printf(msg_verbose,
	     "Assigning density field on grid with template n_mas=%d\n",
	     f.n_mas);

  grid.clear();
  grid.boxsize= boxsize;

  const double dx_inv = grid.nc/boxsize;

  double* const d = grid.fx;

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
    w2_sum += w2;    
    nw2_sum += nbar*w2;
    
    double x[] = {(p->x[0] - x0[0])*dx_inv,
		  (p->x[1] - x0[1])*dx_inv,
		  (p->x[2] - x0[2])*dx_inv};
    f(x, w, grid);
  }

  grid.total_weight = w_sum;
  grid.raw_noise = w2_sum;
  grid.normalisation = nw2_sum;
  grid.n_mas = f.n_mas;
  
  // time duration
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time mass_assignment template %le\n",
	     std::chrono::duration<double>(te - ts).count());
}

#endif
