#include <cmath>
#include <cassert>
#include "config.h"
#include "error.h"
#include "msg.h"
#include "grid.h"
#include "power_spectrum.h"

PowerSpectrum::PowerSpectrum(const double kmin_, const double kmax_,
			     const double dk_) :
  kmin(kmin_), kmax(kmax_), dk(dk_), dk_inv(1.0/dk)
{
  n = (int) ceil((kmax - kmin)/dk);
  nmode_hist= (int*) malloc(sizeof(int)*n); assert(nmode_hist);
  k_hist= (double*) malloc(sizeof(double)*n*4); assert(k_hist);

  P0_hist= k_hist + n;
  P2_hist= P0_hist + n;
  P4_hist= P2_hist + n;
}

PowerSpectrum::~PowerSpectrum()
{
  free(nmode_hist);
  free(k_hist);
}

void PowerSpectrum::clear()
{
  for(int i=0; i<n; ++i) {
    nmode_hist[i]= 0;
    k_hist[i]  = 0.0;
    P0_hist[i] = 0.0;
    P2_hist[i] = 0.0;
    P4_hist[i] = 0.0;
  }
}

void PowerSpectrum::normalise()
{
  for(int i=0; i<n; ++i) {
    if(nmode_hist[i] == 0) continue;

    k_hist[i]  /= nmode_hist[i];
    P0_hist[i] /= nmode_hist[i];
    P2_hist[i] /= nmode_hist[i];
    P4_hist[i] /= nmode_hist[i];
  }
}

static inline float w(const float tx)
{
  return tx == 0.0f ? 1.0f : sin(tx)/tx;
}

void power_spectrum_compute_multipoles(Grid const * const grid,
			    const Float nbar, const Float neff,
			    PowerSpectrum* const ps)
{
  // Power spectrum subtracting shot noise and aliasing correction
  // P(k) ~ k^neff
  const size_t nc= grid->nc;
  const Float boxsize= grid->boxsize;

  if(nc == 0) {
    msg_printf(msg_fatal, "Error: grid->nc is zero.");
    throw AssertionError();
  }

  if(boxsize <= 0.0) {
    msg_printf(msg_fatal, "Error: grid->boxsize is not positive.");
    throw AssertionError();
  }
  //const Float kmin, const Float kmax, const Float dk)
			    
  const Float knq= (M_PI*nc)/grid->boxsize;
  const Float knq_inv= boxsize/(M_PI*nc);
  const Float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const Float sin_fac= 0.5*boxsize/nc;

  const Float fac= 2.0*M_PI/boxsize;
  const size_t ncz= nc/2+1;

  const Float kp= 2.0*M_PI*nc/boxsize;

  const Float nbar_inv= nbar > 0.0f ? 1.0/nbar : 0.0f;

  const int na=2;

  //  if(kmax == 0.0f) kmax= nc*M_PI/boxsize; //Nyquist frequency

  Complex* const delta_k= grid->fk;

  //Histogram Pgg(kmin, kmax, dk);  // 10.0/boxsize
  //Histogram P2(kmin, kmax, dk);  // quadrupole
  //Histogram P4(kmin, kmax, dk);  // hexapole
  
  for(int ix=0; ix<nc; ++ix) {
    Float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    Float sintx= sin(sin_fac*kx);

    float c1x= 1.0f - 2.0f/3.0f*sintx*sintx;

    for(int iy=0; iy<nc; ++iy) {
      Float ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);

      Float sinty= sin(sin_fac*ky);
      Float c1y= 1.0 - 2.0/3.0*sinty*sinty;

      int iz0 = !(kx > 0.0f || (kx == 0.0 && ky > 0.0));

      // Avoid double counting on kz=plain
      // k=(0,0,0) dropped because this is 0
      // iz0= 0 if kx>0 or (kx == 0 and ky > 0)
      //      1 otherwize
      //
      
      for(int iz=iz0; iz<nc/2+1; ++iz) {
	Float kz= fac*iz;
	Float sintz= sin(sin_fac*kz);
	Float c1z= 1.0f - 2.0f/3.0f*sintz*sintz;

	Float k= sqrt(kx*kx + ky*ky + kz*kz);
	Float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005
	Float mu2= kz*kz/(kx*kx + ky*ky + kz*kz);
	
	if(k <= knq) {

	// C2 function in Jing 2005
	Float c2gg= 0.0f, c2gu= 0.0f;
	for(int ax=-na; ax<na; ++ax) {
	  Float kax= kx + ax*kp;
	  for(int ay=-na; ay<na; ++ay) {
	    Float kay= ky + ay*kp;
	    for(int az=-na; az<na; ++az) {
	      Float kaz= kz + az*kp;
	      Float ka= sqrt(kax*kax + kay*kay + kaz*kaz);
	      Float kk= ka/k;

	      float w1= w(sin_fac*(kx + ax*kp))*
		        w(sin_fac*(ky + ay*kp))*
                        w(sin_fac*(kz + az*kp));
	      float w2= w1*w1;
	      float w4= w2*w2;
	      float Pkfac= pow(ka/k, neff);;
	      c2gg += w4*Pkfac;
	    }
	  }
	}

	assert(c2gg > 0.0f);
	
	size_t index= ncz*(nc*ix + iy) + iz;
	Float d_re= delta_k[index][0];
	Float d_im= delta_k[index][1];

	// Legendre polynomial
	// P_l = (2 l + 1)/2*int_0^1 P(k) P_l(mu) dmu
	// P_2 = (3 mu^2 + 1)/2
	// P_4 = (35 mu^4 - 30 mu^2 + 3)/8
	float l_2= 7.5f*mu2 - 2.5f;
	float l_4= 1.125*(35.0f*mu2*mu2 - 30.0f*mu2 + 3.0f);

	ps->add(k,
		(pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg,
		l_2*(pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg,
		l_4*(pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg);	
	
	}
      }
    }
  }
}
