//#include <iostream>
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
#include "discrete_multipole.h"

using namespace std;

void multipole_mas_correction_vector(const int nc, const int n_mas,
				     vector<Float>& mas_correction_array)
{
  mas_correction_array.clear();
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

class DiscreteMultipole {
public:
  explicit DiscreteMultipole(Grid const * const grid_,
			     const vector<Float>& discrete_legendre_coef) :
    grid(grid_), coef(discrete_legendre_coef) { }

  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    // Discrete Legendre polynomial
    // computed in discrete_multipole.cpp
    double l2= coef[5*ik] + coef[5*ik + 1]*mu2;
    double l4= coef[5*ik + 2] + coef[5*ik + 3]*mu2 + coef[5*ik + 4]*mu2*mu2;
    double delta2= norm(grid->fk[index])*corr;

    P.p0[ik] += delta2;
    P.p2[ik] += l2*delta2;
    P.p4[ik] += l4*delta2;
  }
  Grid const * const grid;
  const vector<Float>& coef;
};

// Multipole of P(k) grid, not delta(k) grid
class PowerMultipole {
public:
  explicit PowerMultipole(Grid const * const grid_) :
    grid(grid_) { }

  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    double l2= 7.5*mu2 - 2.5;
    double l4= (1.125*35.0)*mu2*mu2 - (1.125*30.0)*mu2 + (1.125*3.0);
    double pk= real(grid->fk[index])*corr;
    P.p0[ik] += pk;
    P.p2[ik] += l2*pk;
    P.p4[ik] += l4*pk;
  }
  Grid const * const grid;
};

// Discrete Legendre Multipole of P(k) grid, not delta(k) grid
class PowerDiscreteMultipole {
public:
  explicit PowerDiscreteMultipole(Grid const * const grid_,
			  const vector<Float>& discrete_legendre_coef) :
    grid(grid_), coef(discrete_legendre_coef) { }

  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    // Discrete Legendre polynomial
    // computed in discrete_multipole.cpp

    
    double l2= coef[5*ik] + coef[5*ik + 1]*mu2;
    double l4= coef[5*ik + 2] + coef[5*ik + 3]*mu2 + coef[5*ik + 4]*mu2*mu2;
    double pk= real(grid->fk[index])*corr; // grid of P(k)

    P.p0[ik] += pk;
    P.p2[ik] += l2*pk;
    P.p4[ik] += l4*pk;
  }
  Grid const * const grid;
  const vector<Float>& coef;
};


class MultipoleBianchi {
public:
  explicit MultipoleBianchi(Grid const * const grid_,
			    Grid const * const grid0_,
			    Grid const * const grid2_,
			    Grid const * const grid4_) :
    grid(grid_), grid0(grid0_), grid2(grid2_), grid4(grid4_) {
    assert(grid0->nc == grid->nc);
    assert(grid2->nc == grid->nc);
    assert(grid4->nc == grid->nc);
  }
  
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    Complex delta_k = grid->fk[index];
    Complex delta0_k = grid0->fk[index];
    Complex delta2_k = grid2->fk[index] - 0.5*delta0_k;
    Complex delta4_k = grid4->fk[index] - 2.5*delta2_k - 7.0/8.0*delta0_k;


    P.p0[ik] +=     (  delta0_k.real()*delta_k.real()
		     + delta0_k.imag()*delta_k.imag())*corr;
    P.p2[ik] += 5.0*(  delta2_k.real()*delta_k.real() 
	             + delta2_k.imag()*delta_k.imag())*corr;
    P.p4[ik] += 9.0*(  delta4_k.real()*delta_k.real()
	             + delta4_k.imag()*delta_k.imag())*corr;
  }

  Grid const * const grid;
  Grid const * const grid0;
  Grid const * const grid2;
  Grid const * const grid4;
};



class MultipoleScoccimarro {
 public:
  explicit MultipoleScoccimarro(Grid const * const grid_,
				Grid const * const grid0_,
				Grid const * const grid2_) :
    grid(grid_), grid0(grid0_), grid2(grid2_) {
    assert(grid2->nc == grid->nc);
  }
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    Complex delta_k = grid->fk[index];
    Complex delta0_k = grid0->fk[index];
    Complex delta2_k = grid2->fk[index] - 0.5*delta0_k;

    double P0_hat = delta0_k.real()*delta_k.real()
                    + delta0_k.imag()*delta_k.imag()*corr;
    double P2_hat = 5.0*(delta2_k.real()*delta_k.real() 
	            + delta2_k.imag()*delta_k.imag())*corr;

    double delta2_sq = norm(delta2_k)*corr;
    //abort(); //DEBUG!!!!!
	
    P.p0[ik] += P0_hat;
    P.p2[ik] += P2_hat;
    P.p4[ik] += 17.5*delta2_sq - P2_hat - 3.5*P0_hat;
  }
  Grid const * const grid;
  Grid const * const grid0;
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
    
    // Imaginary part of delta1 delta0^*
    // (2l + 1) = 3.0 is multiplied here
    double P1_hat = 3.0*(  delta1_k.imag()*delta0_k.real() 
			 - delta1_k.real()*delta0_k.imag())*corr;

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

    // 2l + 1 = 3 for l = 1
    double P1_hat = 3.0*(  delta1_k.imag()*delta0_k.real()
			 - delta1_k.real()*delta0_k.imag())*corr;

    // 2l + 1 = 2*3 + 1 = 7
    double P3_hat = 7.0*(  delta3_k.imag()*delta0_k.real() 
			 - delta3_k.real()*delta0_k.imag())*corr;

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
multipole_compute_discrete_multipoles(const bool is_delta,
				      const double k_min, const double k_max,
				      const double dk, 
				      Grid const * const grid,
				      const bool subtract_shotnoise,
				      const bool correct_mas,
				      const int line_of_sight)
{
  // is_delta = true:  grid is delta(k)
  //            false: grid is P(k)
  vector<double> coef;
  discrete_multipole_compute_legendre(k_min, k_max, dk, grid->boxsize,
				      coef);

  if(is_delta)
    return compute_multipoles_template(k_min, k_max, dk, 
				       DiscreteMultipole(grid, coef),
				       subtract_shotnoise, correct_mas,
				       line_of_sight);

  return compute_multipoles_template(k_min, k_max, dk, 
				     PowerDiscreteMultipole(grid, coef),
				     subtract_shotnoise, correct_mas,
				     line_of_sight);

}

PowerSpectrum*
multipole_compute_power_multipoles(const double k_min, const double k_max,
				   const double dk, 
				   Grid const * const grid,
				   const bool subtract_shotnoise,
				   const bool correct_mas,
				   const int line_of_sight)
{
  return compute_multipoles_template(k_min, k_max, dk, 
				     PowerMultipole(grid),
				     subtract_shotnoise, correct_mas,
				     line_of_sight);
}

PowerSpectrum*
multipole_compute_yamamoto_scoccimarro(const double k_min, const double k_max,
				       const double dk,
				       Grid const * const grid,
				       Grid const * const grid0,
				       Grid const * const grid2,
				       const bool subtract_shotnoise,
				       const bool correct_mas)
{
  return compute_multipoles_template(k_min, k_max, dk,
				     MultipoleScoccimarro(grid, grid0, grid2),
				     subtract_shotnoise, correct_mas);
}

PowerSpectrum*
multipole_compute_yamamoto_bianchi(const double k_min, const double k_max, const double dk, 
			   Grid const * const grid,
			   Grid const * const grid0,
			   Grid const * const grid2,
			   Grid const * const grid4,
			   const bool subtract_shotnoise,
			   const bool correct_mas)
{
  return compute_multipoles_template(k_min, k_max, dk,
				    MultipoleBianchi(grid, grid0, grid2, grid4),
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

