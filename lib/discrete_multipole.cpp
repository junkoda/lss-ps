#include <vector>
#include <cmath>
#include "config.h"
#include "discrete_multipole.h"

using namespace std;

struct DiscruteLegendre {
  vector<Float> coef;
};

void discrete_multipole_compute_legendre(const Float k_min, const Float k_max, const Float dk, const Float boxsize, vector<Float>& coef)
{
  // Compute coefficients of discrete legendre polynomials
  // k_min   (Float): [h/Mpc]
  // k_max   (Float): [h/Mpc]
  // dk      (Float): [h/Mpc] bin width
  // boxsize (Float): Periodic box length on a side
  //
  // Result
  //   coef[5*ibin + iparam] a0(2) a2(2) a0(4) a2(4) a4(4)
  //
  // Discrete Legendre polynomials
  // P0(x) = 1
  // P2(x) = a_0(2) + a_2(2) x^2
  // P4(x) = a_0(4) + a_2(2) x^2 + a_4(4) x^4

  
  const Float fac= 2*M_PI/boxsize;
  const int ik_max= ceil(k_max/fac);


  const int nbin= ceil((k_max - k_min)/dk);

  vector<int> nmodes(nbin);
  vector<Float> mu2(nbin), mu4(nbin), mu6(nbin), mu8(nbin);
  
  for(int ikx=-ik_max; ikx<=ik_max; ++ikx) {
    for(int iky=-ik_max; iky<=ik_max; ++iky) {
      // skip half of the ikz=0 plane
      int iz0 = !(ikx > 0 || (ikx == 0 && iky > 0));
      for(int ikz=iz0; ikz<=ik_max; ++ikz) {
	Float ik2= static_cast<Float>(ikx*ikx + iky*iky + ikz*ikz);
	Float m2= ikz/ik2;
	Float ik= sqrt(ik2);

	int ibin= floor((fac*ik - k_min)/dk);
	if(0 <= ibin && ibin < nbin) {
	  nmodes[ibin]++;
	  mu2[ibin] += m2;
	  mu4[ibin] += m2*m2;
	  mu6[ibin] += m2*m2*m2;
	  mu8[ibin] += m2*m2*m2*m2;
	}
      }
    }
  }

  coef.resize(5*nbin, 0.0);
  
  for(int ibin=0; ibin<nbin; ++ibin) {
    if(nmodes[ibin] > 0) {
      mu2[ibin] /= nmodes[ibin];
      mu4[ibin] /= nmodes[ibin];
      mu6[ibin] /= nmodes[ibin];
      mu8[ibin] /= nmodes[ibin];
      
      Float a2= 1.0/mu2[ibin];
      Float norm= sqrt(5.0*(1.0 + 2.0*a2*mu2[ibin] + mu4[ibin]));

      coef[5*ibin    ]= -1.0/norm;    // a_0(2)
      coef[5*ibin + 1]=  a2/norm;      // a_2(2)

      Float det= mu2[ibin]*mu6[ibin] - mu4[ibin]*mu4[ibin];
      a2= (mu6[ibin] - mu4[ibin])/det;
      Float a4= (mu2[ibin] - mu4[ibin])/det;
      norm = sqrt(1.0 + a2*a2*mu4[ibin] + a4*a4*mu8[ibin]
	          + 2.0*(a2*mu2[ibin] + a4*mu4[ibin] + a2*a4*mu6[ibin]));
      Float sign= a4 > 0 ? 1.0 : -1.0;
      coef[5*ibin + 2]= sign/norm;     // a_0(4)
      coef[5*ibin + 3]= sign*a2/norm;  // a_2(4)
      coef[5*ibin + 4]= sign*a4/norm;  // a_4(4)
    }
  }

}
