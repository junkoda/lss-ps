#include <iostream>
#include <chrono>
#include <cmath>
#include <cassert>
#include "config.h"
#include "error.h"
#include "msg.h"
#include "grid.h"
#include "power_spectrum.h"

//
// class Power Spectrum
//
PowerSpectrum::PowerSpectrum(const double kmin_, const double kmax_,
			     const double dk_) :
  kmin(kmin_), kmax(kmax_), dk(dk_), dk_inv(1.0/dk)
{
  n = (int) round((kmax - kmin)/dk);
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

inline float w(const float tx)
{
  return tx == 0.0f ? 1.0f : sin(tx)/tx;
}

//
// Compute power spectrum from Fourier grid aka shell average
//

void power_spectrum_compute_multipoles(Grid const * const grid,
				       PowerSpectrum* const ps,
				       const bool subtract_shotnoise)
{
  auto ts = std::chrono::high_resolution_clock::now();
  
  // Power spectrum subtracting shot noise and aliasing correction
  // P(k) ~ k^neff (unused now)
  const size_t nc= grid->nc;
  const Float boxsize= grid->boxsize;
  ps->clear();

  msg_printf(msg_verbose,
	     "computing power spectrum multipoles.\n");
  
  if(nc == 0) {
    msg_printf(msg_error, "Error: grid->nc is zero.\n");
    throw AssertionError();
  }

  if(boxsize <= 0.0) {
    msg_printf(msg_error, "Error: grid->boxsize is not positive.\n");
    throw AssertionError();
  }

  if(grid->mode != grid_fourier_space) {
    msg_printf(msg_error, "Error: grid is not in Fourier space.\n");
    throw AssertionError();
  }
			    
  const Float knq= (M_PI*nc)/grid->boxsize;
  const Float knq_inv= boxsize/(M_PI*nc);
  //const Float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  std::cerr << "Total weight " << grid->total_weight << std::endl;
  //const Float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 3)/grid->total_weight;
  const Float pk_fac= boxsize*boxsize*boxsize/(grid->total_weight*grid->total_weight);
  const Float sin_fac= 0.5*boxsize/nc;

  const Float fac= 2.0*M_PI/boxsize;
  const size_t nckz= nc/2+1;

  const Float kp= 2.0*M_PI*nc/boxsize;

  const Float nbar_inv= subtract_shotnoise ? grid->shot_noise : 0;
  msg_printf(msg_info, "Shot noise subtraction: %e\n", nbar_inv);

  Complex* const delta_k= grid->fk;

  for(size_t ix=0; ix<nc; ++ix) {
    Float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    //Float sintx= sin(sin_fac*kx);

    // ToDo: this is for no interlacing derive a correct one
    //float c1x= 1.0 - 2.0/3.0*sintx*sintx;

    for(size_t iy=0; iy<nc; ++iy) {
      Float ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);

      //Float sinty= sin(sin_fac*ky);
      //Float c1y= 1.0 - 2.0/3.0*sinty*sinty;

      size_t iz0 = !(kx > 0.0f || (kx == 0.0 && ky > 0.0));

      // Avoid double counting on kz=plain
      // k=(0,0,0) dropped because this is 0
      // iz0= 0 if kx>0 or (kx == 0 and ky > 0)
      //      1 otherwize
      //
      
      for(size_t iz=iz0; iz<nckz; ++iz) {
	Float kz= iz <= nc/2 ? fac*iz : fac*(iz-nc);
	//Float sintz= sin(sin_fac*kz);
	//Float c1z= 1.0 - 2.0/3.0*sintz*sintz;

	Float k= sqrt(kx*kx + ky*ky + kz*kz);
	// This shot_noise is without interlacing
	//Float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005
	Float shot_noise= nbar_inv;
	Float mu2= kz*kz/(kx*kx + ky*ky + kz*kz);
	
	if(k <= knq) {

	  Float w1= w(sin_fac*(kx))*
	            w(sin_fac*(ky))*
	            w(sin_fac*(kz));
	  Float w2= w1*w1;
	  Float c2gg= w2*w2;

	
	  size_t index= nckz*(nc*ix + iy) + iz;
	  Float d_re= delta_k[index][0];
	  Float d_im= delta_k[index][1];

	// Legendre polynomial
	// P_l = (2 l + 1)/2*int_0^1 P(k) P_l(mu) dmu
	// P_2 = (3 mu^2 + 1)/2
	// P_4 = (35 mu^4 - 30 mu^2 + 3)/8
	Float l_2= 7.5f*mu2 - 2.5f;
	Float l_4= 1.125*(35.0f*mu2*mu2 - 30.0f*mu2 + 3.0f);
	Float delta2= pk_fac*(d_re*d_re + d_im*d_im);

	//if(ix < nc/2 && iy < nc/2 && iz < nc/2) // debug!!!
	ps->add(k,
		delta2/c2gg, 
		l_2*(delta2/c2gg - shot_noise),
		l_4*(delta2/c2gg - shot_noise));	
	
	}
      }
    }
  }
  ps->normalise();

  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time multipoles %le\n",
	     std::chrono::duration<double>(te - ts).count());

}

void power_spectrum_compute_multipoles_jing(Grid const * const grid,
					    PowerSpectrum* const ps,
					    const bool subtract_shotnoise,
					    const Float neff)

{
  // Power spectrum subtracting shot noise (optional) and aliasing correction
  // using Jing (2005) method
  // P(k) ~ k^neff
  auto ts = std::chrono::high_resolution_clock::now();
  
  const size_t nc= grid->nc;
  const Float boxsize= grid->boxsize;
  ps->clear();

  msg_printf(msg_verbose, "computing power spectrum multipoles.\n");
  msg_printf(msg_info, "Jing alias correction.\n");
  
  if(nc == 0) {
    msg_printf(msg_error, "Error: grid->nc is zero.\n");
    throw AssertionError();
  }

  if(boxsize <= 0.0) {
    msg_printf(msg_error, "Error: grid->boxsize is not positive.\n");
    throw AssertionError();
  }

  if(grid->mode != grid_fourier_space) {
    msg_printf(msg_error, "Error: grid is not in Fourier space.\n");
    throw AssertionError();
  }
			    
  const Float knq= (M_PI*nc)/grid->boxsize;
  const Float knq_inv= boxsize/(M_PI*nc);
  const Float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const Float sin_fac= 0.5*boxsize/nc;

  const Float fac= 2.0*M_PI/boxsize;
  const size_t ncz= nc/2+1;

  const Float kp= 2.0*M_PI*nc/boxsize;

  const Float nbar_inv= subtract_shotnoise ? grid->shot_noise : 0;
  msg_printf(msg_info, "Shot noise subtraction: %e\n", nbar_inv);
  
  const int na=2;

  Complex* const delta_k= grid->fk;

  for(int ix=0; ix<nc; ++ix) {
    Float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    Float sintx= sin(sin_fac*kx);

    float c1x= 1.0 - 2.0/3.0*sintx*sintx;

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
	Float c1z= 1.0 - 2.0/3.0*sintz*sintz;

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
  ps->normalise();

  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time multipoles Jing %le\n",
	     std::chrono::duration<double>(te - ts).count());

}

/*
void power_spectrum_compute_shotnoise(Grid const * const grid,
				      PowerSpectrum* const ps)

{
  // Correct power spectrum aliasing with Jing C1 function
  // which assumes that the power spectrum is flat (white noise)
  const size_t nc= grid->nc;
  const Float boxsize= grid->boxsize;
  ps->clear();

  msg_printf(msg_verbose, "computing shotnoise power spectrum.\n");
  
  if(nc == 0) {
    msg_printf(msg_error, "Error: grid->nc is zero.\n");
    throw AssertionError();
  }

  if(boxsize <= 0.0) {
    msg_printf(msg_error, "Error: grid->boxsize is not positive.\n");
    throw AssertionError();
  }

  if(grid->mode != grid_fourier_space) {
    msg_printf(msg_error, "Error: grid is not in Fourier space.\n");
    throw AssertionError();
  }
			    
  const Float knq= (M_PI*nc)/grid->boxsize;
  const Float knq_inv= boxsize/(M_PI*nc);
  const Float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const Float sin_fac= 0.5*boxsize/nc;

  const Float fac= 2.0*M_PI/boxsize;
  const size_t ncz= nc/2+1;

  const Float kp= 2.0*M_PI*nc/boxsize;

  Complex* const delta_k= grid->fk;

  for(int ix=0; ix<nc; ++ix) {
    Float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    Float sintx= sin(sin_fac*kx);

    float c1x= 1.0 - 2.0/3.0*sintx*sintx;

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
	Float c1z= 1.0 - 2.0/3.0*sintz*sintz;

	Float k= sqrt(kx*kx + ky*ky + kz*kz);
	Float shot_noise= c1x*c1y*c1z; // C1 function in Jing 2005
	
	if(k <= knq) {

	size_t index= ncz*(nc*ix + iy) + iz;
	Float d_re= delta_k[index][0];
	Float d_im= delta_k[index][1];

	ps->add(k,
		pk_fac*(d_re*d_re + d_im*d_im)/shot_noise, 0.0, 0.0);	
	}
      }
    }
  }
  ps->normalise();
}

void power_spectrum_compute_shotnoise(const size_t nc, const Float boxsize,
				      const Float nbar,
				      PowerSpectrum* const ps)
{
  // Power spectrum without any correction
  ps->clear();

  msg_printf(msg_verbose, "computing shotnoise multipoles.\n");
  msg_printf(msg_info, "No alias correction\n");
  
  if(nc == 0) {
    msg_printf(msg_error, "Error: nc is zero.\n");
    throw AssertionError();
  }

  if(boxsize <= 0.0) {
    msg_printf(msg_error, "Error: boxsize is not positive.\n");
    throw AssertionError();
  }

  if(nbar == 0.0) {
    msg_printf(msg_error, "Error: nbar must be positive.\n");
    throw AssertionError();
  }

  const Float knq= (M_PI*nc)/boxsize;
  const Float knq_inv= boxsize/(M_PI*nc);
  const Float sin_fac= 0.5*boxsize/nc;

  const Float fac= 2.0*M_PI/boxsize;
  const size_t ncz= nc/2+1;

  const Float kp= 2.0*M_PI*nc/boxsize;

  const Float nbar_inv= 1.0/nbar;

  for(int ix=0; ix<nc; ++ix) {
    Float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    Float sintx= sin(sin_fac*kx);
    
    float c1x= 1.0 - 2.0/3.0*sintx*sintx;

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
	Float c1z= 1.0 - 2.0/3.0*sintz*sintz;

	Float w1= w(sin_fac*(kx))*
	          w(sin_fac*(ky))*
	          w(sin_fac*(kz));
	Float w2= w1*w1;
	Float c2gg= w2*w2;

	Float k= sqrt(kx*kx + ky*ky + kz*kz);
	Float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005
	Float mu2= kz*kz/(kx*kx + ky*ky + kz*kz);
	
	if(k <= knq) {

	// Legendre polynomial
	// P_l = (2 l + 1)/2*int_0^1 P(k) P_l(mu) dmu
	// P_2 = (3 mu^2 + 1)/2
	// P_4 = (35 mu^4 - 30 mu^2 + 3)/8
	float l_2= 7.5f*mu2 - 2.5f;
	float l_4= 1.125*(35.0f*mu2*mu2 - 30.0f*mu2 + 3.0f);

	ps->add(k,
		shot_noise,
		l_2*shot_noise,
		l_4*shot_noise);	
	
	}
      }
    }
  }
  ps->normalise();
}
*/
