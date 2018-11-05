#include "yamamoto.h"

//
// Yamamoto-Bianchi function object for compute_multipoles() algorithm
//
class MultipoleBianchi {
 public:
  explicit MultipoleBianchi(GridFFT const * const grid_delta0,
			    GridFFT const * const grid_delta2,
			    GridFFT const * const grid_delta4) :
    delta0(grid_delta0),
    delta2(grid_delta2),
    delta4(grid_delta4),

    assert(delta2->nc == delta0->nc);
    assert(delta4->nc == delta0->nc);
  }
  
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    std::complex<double> delta0_k= delta0->fk[index];

    std::complex<double> delta2_k= delta2->fk[index] - 0.5*delta0_k;
    
    std::complex<double> delta4_k= delta4->fk[index] - 2.5*delta2_k
                                   - 7.0/8.0*delta0_k;


    P.p0[ik] += norm(delta0_k)*corr;
    P.p2[ik] += 5.0*(delta2_k.real()*delta0_k.real() 
	           + delta2_k.imag()*delta0_k.imag())*corr;
    P.p4[ik] += 9.0*(delta4_k.real()*delta0_k.real()
	           + delta4_k.imag()*delta0_k.imag())*corr;
  }

  GridFFT const * const delta0;
  GridFFT const * const delta2;
  GridFFT const * const delta4;	   
};


//
// Yamamoto-Scoccimarro function object for compute_multipoles() algorithm
//
class MultipoleScoccimarro {
 public:
  MultipoleScoccimarro(Grid const * const delta0_,
		       Grid const * const delta2_) :
    delta0(delta0_), delta2(delta2_) {
    assert(delta0->nc == delta2->nc);
  }
  void operator()(const size_t index, const double mu2, const double corr,
		  const int ik, PowerSpectrum& P) const {
    std::complex<double> delta0_k = delta0->fk[index];
    std::complex<double> delta2_k = delta2->fk[index] - 0.5*delta0_k;
    
    double P0_hat = norm(delta0_k)*corr;
    double P2_hat = 5.0*(delta2_k.real()*delta0_k.real() 
	               + delta2_k.imag()*delta0_k.imag())*corr;

    double delta2_sq = norm(delta2_k)*corr;
	
    P.p0[ik] += P0_hat;
    P.p2[ik] += P2_hat;
    P.p4[ik] += 17.5*delta2_sq - P2_hat - 3.5*P0_hat;
  }
  Grid const * const delta0;
  Grid const * const delta2;
};


