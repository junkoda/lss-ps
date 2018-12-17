#ifndef MULTIPOLE_H
#define MULTIPOLE_H 1

#include <vector>
#include <algorithm>
#include <chrono>

#include "power_spectrum.h"
#include "grid.h"
#include "mas_correction.h"

PowerSpectrum*
multipole_compute_monopole(const double k_min, const double k_max,
			   const double dk,
			   Grid const * const grid,
			   const bool subtract_shotnoise,
			   const bool correct_mas);

PowerSpectrum*
multipole_compute_plane_parallel(const double k_min, const double k_max,
				 const double dk, 
				 Grid const * const grid,
				 const bool subtract_shotnoise,
				 const bool correct_mas,
				 const int line_of_sight);

PowerSpectrum*
multipole_compute_discrete_multipoles(const bool is_delta,
				      const double k_min, const double k_max,
				      const double dk, 
				      Grid const * const grid,
				      const bool subtract_shotnoise,
				      const bool correct_mas,
				      const int line_of_sight);

PowerSpectrum*
multipole_compute_power_multipoles(const double k_min, const double k_max,
				   const double dk, 
				   Grid const * const grid,
				   const bool subtract_shotnoise,
				   const bool correct_mas,
				   const int line_of_sight);

PowerSpectrum*
multipole_compute_yamamoto_scoccimarro(const double k_min, const double k_max,
			   const double dk,
			   Grid const * const grid,
			   Grid const * const grid0,
			   Grid const * const grid2,
			   const bool subtract_shotnoise,
			   const bool correct_mas);

PowerSpectrum*
multipole_compute_yamamoto_bianchi(const double k_min, const double k_max,
			   const double dk,
			   Grid const * const grid,
			   Grid const * const grid0,
			   Grid const * const grid2,
			   Grid const * const grid4,
			   const bool subtract_shotnoise,
			   const bool correct_mas);

PowerSpectrum*
multipole_compute_yamamoto_odd_multipoles(const double k_min,
					  const double k_max, const double dk, 
					  Grid const * const grid,
					  Grid const * const grid1,
					  Grid const * const grid3,
					  const bool subtract_shotnoise,
					  const bool correct_mas);

PowerSpectrum*
multipole_compute_yamamoto_zero(const double k_min,
				const double k_max, const double dk, 
				Grid const * const grid,
				Grid const * const grid0,
				const bool subtract_shotnoise,
				const bool correct_mas);


//
// The main algorithm for multipole
//
template <typename F>
PowerSpectrum* compute_multipoles_template(const double k_min,
					   const double k_max,
					   const double dk, 
					   const F f,
					   const bool subtract_shotnoise,
					   const bool correct_mas,
					   const int los=2)
{
  auto ts = std::chrono::high_resolution_clock::now();

  assert(0 <= los && los < 3);

  Grid const * const grid = f.grid; assert(grid);
  const int nc= grid->nc;
  const double boxsize= grid->boxsize;
  const int n_mas= correct_mas ? grid->n_mas : 0;

  std::vector<Float>& mas_correction_array=
    mas_correction_get_array(nc, n_mas);

  assert(nc > 0);
  assert(boxsize > 0.0);
  assert(grid->mode == GridMode::fourier_space);
  assert(mas_correction_array.size() == (size_t) nc);
  

  // binning
  const double k_fundamental= 2.0*M_PI/boxsize;
  const double k_nq= M_PI/boxsize*nc;
  const int n= (int) round((k_max - k_min)/dk);
  const double ik_min= k_min/k_fundamental;
  const double ik_max= std::min(k_nq, k_max)/k_fundamental;
  const double idk= dk/k_fundamental;

  PowerSpectrum* const ps= new PowerSpectrum(n);
  PowerSpectrum& P= *ps;

  const size_t nckz= nc/2+1;
  const int ik_nq= nc/2;
  const int ik_max2 = ik_max*ik_max;

  const double shot_noise= subtract_shotnoise ? grid->shot_noise : 0;
  msg_printf(msg_info, "Shot noise subtraction: %e\n", shot_noise);
  msg_printf(msg_info, "Pk normalisation %e\n", grid->pk_normalisation);

#ifdef _OPENMP
  #pragma omp parallel num_threads(omp_get_max_threads())
#endif
  {
    PowerSpectrum ps_local(n);

#ifdef _OPENMP
    #pragma omp for
#endif
    for(int ix=0; ix<nc; ++ix) {
      int k[3];
      k[0] = ix <= ik_nq ? ix : ix - nc;
      if(k[0] >= ik_max || k[0] <= -ik_max)
	continue; 
      
      double corr_x = mas_correction_array[ix];
      
      for(int iy=0; iy<nc; ++iy) {
	k[1] = iy <= ik_nq ? iy : iy - nc;
	int kk = k[0]*k[0] + k[1]*k[1];
	if(kk > ik_max2) continue;
	
	double corr_xy = corr_x * mas_correction_array[iy];
	
	int kz0 = !(k[0] > 0 || (k[0] == 0 && k[1] > 0));
	
	// Avoid double counting on kz=0 plain
	// k=(0,0,0) dropped because this is 0
	// iz0= 0 if kx>0 or (kx == 0 and ky > 0)
	//      1 otherwize
	//
	
	for(int iz=kz0; iz<ik_nq; ++iz) {
	  k[2]= iz;
	  double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
	  if(k2 > ik_max2) break;
	  
	  double kmag= sqrt(k2);
	  int i= static_cast<int>(floor((kmag - ik_min)/idk));
	
	  if(0 <= i && i < n) {
	    size_t index= nckz*(nc*ix + iy) + iz;
	    double mu2= (k[los]*k[los])/k2;	
	    double corr_xyz = corr_xy * mas_correction_array[iz];

	    ps_local.nmodes[i]++;
	    ps_local.k[i] += kmag; 

	    f(index, mu2, corr_xyz, i, ps_local);
	  }	
	}
      }
    }
			   

    /// Local power spectrum in each core is added up to total
#ifdef _OPENMP
    #pragma omp critical (__COLLECT_PK_1D__)
#endif
    {
      P += ps_local;
    }
  }			   

  // Normalize by the number of modes
  const double pk_fac= grid->pk_normalisation; 
  for(int i=0; i<n; ++i) {
    if(P.nmodes[i] > 0) {
      P.k[i] *= k_fundamental/P.nmodes[i];
      
      double fac = pk_fac/P.nmodes[i];
      P.p0[i] = fac*P.p0[i] - shot_noise;
      P.p1[i] = fac*P.p1[i];
      P.p2[i] = fac*P.p2[i];
      P.p3[i] = fac*P.p3[i];
      P.p4[i] = fac*P.p4[i];
    }
  }
  
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_info, "Time multipoles %e\n",
             std::chrono::duration<double>(te - ts).count());

  ps->shot_noise= shot_noise;

  return ps;
}


#endif
