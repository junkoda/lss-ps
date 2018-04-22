#include <iostream>
#include <vector>
#include <complex>
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
 /*
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
 */


vector<Float> mass_assignment_correction_array(const int nc, const int deg)
{
  vector<Float> v(nc, 1.0);
  
  if(deg == 0) // mas_correction = false
    return v;
  
  const int knq = nc/2;
  const Float fac= M_PI/nc;

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i=1; i<nc; ++i) {
    int k= i <= knq ? i : i - nc;
    Float sinc = sin(fac*k)/(fac*k);
    v[i] = 1.0/pow(sinc, deg);
  }
  
  return v;
}

void mass_assignment_correct_mas(Grid* const grid)
{
  auto ts = std::chrono::high_resolution_clock::now();

  const int nc= grid->nc;
  const double boxsize= grid->boxsize;
  const int n_mas= grid->n_mas;

  assert(nc > 0);
  assert(boxsize > 0.0);
  assert(grid->mode == GridMode::fourier_space);

  vector<Float> mas_correction_array =
    mass_assignment_correction_array(nc, n_mas);

  const size_t nckz= nc/2+1;
  const int ik_nq= nc/2;

  complex<Float>* const fk= grid->fk;

#ifdef _OPENMP
  #pragma omp parallel num_threads(omp_get_max_threads())
#endif
  {
#ifdef _OPENMP
    #pragma omp for
#endif
    for(int ix=0; ix<nc; ++ix) {
      int k[3];
      k[0] = ix <= ik_nq ? ix : ix - nc;
      
      Float corr_x = mas_correction_array[ix];
      
      for(int iy=0; iy<nc; ++iy) {
	k[1] = iy <= ik_nq ? iy : iy - nc;
	
	Float corr_xy = corr_x * mas_correction_array[iy];
	
	int kz0 = !(k[0] > 0 || (k[0] == 0 && k[1] > 0));
	
	for(int iz=kz0; iz<ik_nq; ++iz) {
	  k[2]= iz;
	  
	  Float corr_xyz = corr_xy * mas_correction_array[iz];

	  const size_t index= (ix + iy*nc)*nckz + iz;

	  fk[index] *= corr_xyz;
	  //fk[index] = corr_xyz*fk[index];
	  //fk[index].real() = fk[index].real()/corr_xyz;
	  //fk[index].imag() /= corr_xyz;


	}
      }
    }
  }
			     
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_info, "Time mass assignment correction %e\n",
             std::chrono::duration<double>(te - ts).count());
}



