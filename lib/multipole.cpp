#include <iostream>
#include <vector>
#include <algorithm> // min()
#include <complex>
#include <chrono>
#include <cmath>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "config.h"
#include "msg.h"
#include "multipole.h"
#include "power_spectrum.h"

using namespace std;

// correction factor for mass assignment window function
static int mas_correction_nc= 0;
static int mas_correction_n= -1;
static vector<Float> mas_correction_array;

static void mas_correction_init(const int nc, const int n_mas);

void mas_correction_init(const int nc, const int n_mas)
{
  // initialize mas_correction_array
  if(nc == mas_correction_nc && n_mas == mas_correction_n)
    return;

  // Initialise mas_correction_array
  mas_correction_nc= nc;
  mas_correction_n= nc;
  
  mas_correction_array.clear();
  mas_correction_array.reserve(nc);
  mas_correction_array.assign(nc, (Float) 1);

  if(n_mas == 0) // mas_correction = false
    return;
  
  const int knq = nc/2;
  const Float fac= M_PI/nc;

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i=1; i<nc; ++i) {
    int k= i <= knq ? i : i - nc;
    Float sinc = sin(fac*k)/(fac*k);
    mas_correction_array[i] = 1.0/pow(sinc, 2*n_mas);
  }

  msg_printf(msg_info, "MAS correction array initialised\n");
  
  return;
}


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


  mas_correction_init(nc, n_mas);

  assert(nc > 0);
  assert(boxsize > 0.0);
  assert(grid->mode == GridMode::fourier_space);
  assert(mas_correction_array.size() == (size_t) nc);


  // binning
  const double k_fundamental= 2.0*M_PI/boxsize;
  const double k_nq= M_PI/boxsize*nc;
  const int n= (int) round((k_max - k_min)/dk);
  const double ik_min= k_min/k_fundamental;
  const double ik_max= min(k_nq, k_max)/k_fundamental;
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
      P.p2[i] = fac*P.p2[i];
      P.p4[i] = fac*P.p4[i];
    }
  }
  
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_info, "Time multipoles %e\n",
             std::chrono::duration<double>(te - ts).count());

  ps->shot_noise= shot_noise;

  return ps;
}

//
// Classes for compute_multipoles algorithm
//

class Monopole {
public:
  explicit Monopole(Grid const * const grid_) :
    grid(grid_) { }

  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    P.p0[ik] += norm(grid->fk[index])*corr;
  }
  Grid const * const grid;
};


class Multipole {
public:
  explicit Multipole(Grid const * const grid_) : grid(grid_) { }

  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    // Legendre polynomial
    // (2l + 1) P_l = (2 l + 1)/2*int_0^1 P(k) P_l(mu) dmu
    // (2l + 1) P_2 = 5*(3 mu^2 - 1)/2
    // (2l + 1) P_4 = 9*(35 mu^4 - 30 mu^2 + 3)/8
    double l2= 7.5*mu2 - 2.5;
    double l4= (1.125*35.0)*mu2*mu2 - (1.125*30.0)*mu2 + (1.125*3.0);
    double delta2= norm(grid->fk[index])*corr;
    P.p0[ik] += delta2;
    P.p2[ik] += l2*delta2;
    P.p4[ik] += l4*delta2;
  }
  Grid const * const grid;
};


class MultipoleBianchi {
public:
  explicit MultipoleBianchi(Grid const * const grid_,
			    Grid const * const grid2_,
			    Grid const * const grid4_) :
    grid(grid_), grid2(grid2_), grid4(grid4_) {
    assert(grid2->nc == grid->nc);
    assert(grid4->nc == grid->nc);
  }
  
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    Complex delta0_k = grid->fk[index];
    Complex delta2_k = grid2->fk[index] - 0.5*delta0_k;
    Complex delta4_k = grid4->fk[index] - 2.5*delta2_k - 7.0/8.0*delta0_k;


    P.p0[ik] += norm(delta0_k)*corr;
    P.p2[ik] += 5.0*(delta2_k.real()*delta0_k.real() 
	       + delta2_k.imag()*delta0_k.imag())*corr;
    P.p4[ik] += 9.0*(delta4_k.real()*delta0_k.real()
	       + delta4_k.imag()*delta0_k.imag())*corr;
  }

  Grid const * const grid;
  Grid const * const grid2;
  Grid const * const grid4;
};



class MultipoleScoccimarro {
 public:
  explicit MultipoleScoccimarro(Grid const * const grid_,
				Grid const * const grid2_) :
    grid(grid_), grid2(grid2_) {
    assert(grid2->nc == grid->nc);
  }
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    Complex delta0_k = grid->fk[index];
    Complex delta2_k = grid2->fk[index] - 0.5*delta0_k;

    double P0_hat = norm(delta0_k)*corr;
    double P2_hat = 5.0*(delta2_k.real()*delta0_k.real() 
	            + delta2_k.imag()*delta0_k.imag())*corr;

    double delta2_sq = norm(delta2_k)*corr;
	
    P.p0[ik] += P0_hat;
    P.p2[ik] += P2_hat;
    P.p4[ik] += 17.5*delta2_sq - P2_hat - 3.5*P0_hat;
  }
  Grid const * const grid;
  Grid const * const grid2;
};

//
// Multipole1 (dipole) function object for compute_multipoles() algorithm
//
class Multipole1 {
 public:
  Multipole1(Grid const * const grid0_,
	     Grid const * const grid1_) :
    grid(grid0_), grid1(grid1_) {
  }
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    std::complex<double> delta0_k = grid->fk[index];
    std::complex<double> delta1_k = grid1->fk[index];
    
    double P1_hat = 3.0*(  delta1_k.real()*delta0_k.real() 
		         + delta1_k.imag()*delta0_k.imag())*corr;

    P.p1[ik] += P1_hat;
  }
  Grid const * const grid;
  Grid const * const grid1;
};

//
// Multipole3 (dipole and tripole) function object for
// compute_multipoles() algorithm
//
class Multipole3 {
 public:
  Multipole3(Grid const * const grid0_,
	     Grid const * const grid1_,
	     Grid const * const grid3_) :
    grid(grid0_), grid1(grid1_), grid3(grid3_) {
  }
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    std::complex<double> delta0_k = grid->fk[index];
    std::complex<double> delta1_k = grid1->fk[index];
    std::complex<double> delta3_k = grid3->fk[index] - 1.5*grid1->fk[index] ;

    double P1_hat = 3.0*(  delta1_k.real()*delta0_k.real() 
			 + delta1_k.imag()*delta0_k.imag())*corr;

    double P3_hat = 3.0*(  delta3_k.real()*delta0_k.real() 
		         + delta3_k.imag()*delta0_k.imag())*corr;

    P.p1[ik] += P1_hat;
    P.p3[ik] += P3_hat;
  }
  Grid const * const grid;
  Grid const * const grid1;
  Grid const * const grid3;
};




// Wrappter functions

PowerSpectrum*
multipole_compute_monopole(const double k_min, const double k_max,
			   const double dk, 
			   Grid const * const grid,
			   const bool subtract_shotnoise,
			   const bool correct_mas)
{
  return compute_multipoles_template(k_min, k_max, dk, 
				     Monopole(grid),
				     subtract_shotnoise, correct_mas);
}

PowerSpectrum*
multipole_compute_plane_parallel(const double k_min, const double k_max,
				 const double dk, 
				 Grid const * const grid,
				 const bool subtract_shotnoise,
				 const bool correct_mas,
				 const int line_of_sight)
{
  return compute_multipoles_template(k_min, k_max, dk, 
  				     Multipole(grid),
				     subtract_shotnoise, correct_mas,
				     line_of_sight);
}

PowerSpectrum*
multipole_compute_yamamoto_scoccimarro(const double k_min, const double k_max,
				       const double dk,
				       Grid const * const grid,
				       Grid const * const grid2,
				       const bool subtract_shotnoise,
				       const bool correct_mas)
{
  return compute_multipoles_template(k_min, k_max, dk,
				     MultipoleScoccimarro(grid, grid2),
				     subtract_shotnoise, correct_mas);
}

PowerSpectrum*
multipole_compute_yamamoto_bianchi(const double k_min, const double k_max, const double dk, 
			   Grid const * const grid,
			   Grid const * const grid2,
			   Grid const * const grid4,
			   const bool subtract_shotnoise,
			   const bool correct_mas)
{
  return compute_multipoles_template(k_min, k_max, dk,
				     MultipoleBianchi(grid, grid2, grid4),
				     subtract_shotnoise, correct_mas);
}

PowerSpectrum*
multipole_compute_yamamoto_odd_multipoles(const double k_min,
					  const double k_max, const double dk, 
					  Grid const * const grid,
					  Grid const * const grid1,
					  Grid const * const grid3,
					  const bool subtract_shotnoise,
					  const bool correct_mas)
{
  if(grid3 == nullptr) {
    return compute_multipoles_template(k_min, k_max, dk,
				       Multipole1(grid, grid1),
				       subtract_shotnoise, correct_mas);
  }

  return compute_multipoles_template(k_min, k_max, dk,
				     Multipole3(grid, grid1, grid3),
				     subtract_shotnoise, correct_mas);
}

