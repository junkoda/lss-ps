#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cassert>
#include "config.h"
#include "error.h"
#include "msg.h"
#include "grid.h"
#include "power_spectrum.h"

using namespace std;

//
//
//
static std::vector<Float> window_corr; // array of window function corrections
static int window_corr_n_mas= 0; // n_mas in window_corr
//static void init_window_corr();
void init_window_corr(const size_t nc, const int n_mas);


//
// class Power Spectrum
//
PowerSpectrum::PowerSpectrum(const double kmin_, const double kmax_,
			     const double dk_,
			     const double boxsize) :
  kmin(kmin_), kmax(kmax_), dk(dk_),
  k_fundamental(2.0*M_PI/boxsize),
  ik_min(kmin_/k_fundamental),
  ik_max(kmax_/k_fundamental),
  idk_inv(k_fundamental/dk_),
  n(round((kmax_ - kmin_)/dk_))
{  
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
// 
//
void init_window_corr(const size_t nc, const int n_mas)
{
  if(window_corr.size() == nc && window_corr_n_mas == n_mas)
    return;

  msg_printf(msg_verbose, "Initialising window function correction array.\n");
  window_corr.clear();
  window_corr.reserve(nc);

  if(n_mas == 0) {
    for(int i=0; i<nc; ++i)
      window_corr.push_back(1.0);
    return;
  }

  
  const int knq = nc/2;
  const Float fac= M_PI/nc;
  window_corr.push_back(1.0);
  
  
  for(int i=1; i<nc; ++i) {
    int k= i <= knq ? i : i - nc;
    Float sinc = sin(fac*k)/(fac*k);
    window_corr.push_back(1.0/pow(sinc, 2*n_mas));
  }

  window_corr_n_mas = n_mas;
}

//
// Compute power spectrum from Fourier grid, aka shell average.
//

void power_spectrum_compute_multipoles(Grid const * const grid,
				       PowerSpectrum* const ps,
				       const bool subtract_shotnoise,
				       const bool correct_mas)
{
  auto ts = std::chrono::high_resolution_clock::now();
  
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

  if(correct_mas)
    init_window_corr(nc, grid->n_mas);
  else
    init_window_corr(nc, 0);
  assert(window_corr.size() == nc);
  
  const Float pk_fac= grid->pk_normalisation; 
  cerr << "pk_fac= " << pk_fac << endl;
  const size_t nckz= nc/2+1;
  const int iknq= nc/2;

  const Float shot_noise= subtract_shotnoise ? grid->shot_noise : 0;
  msg_printf(msg_info, "Shot noise subtraction: %e\n", shot_noise);

  Complex* const delta_k= grid->fk;

  for(int ix=0; ix<nc; ++ix) {
    int kx = ix <= iknq ? ix : ix - nc;
    Float corr_x = window_corr[ix];

    for(int iy=0; iy<nc; ++iy) {
      int ky = iy <= iknq ? iy : iy - nc;
      Float corr_xy = corr_x * window_corr[iy];

      int kz0 = !(kx > 0 || (kx == 0 && ky > 0));

      // Avoid double counting on kz=plain
      // k=(0,0,0) dropped because this is 0
      // iz0= 0 if kx>0 or (kx == 0 and ky > 0)
      //      1 otherwize
      //
      
      for(int kz=kz0; kz<iknq; ++kz) {
	Float corr_xyz = corr_xy * window_corr[kz];
	Float k2 = static_cast<double>(kx*kx + ky*ky + kz*kz);
	Float k = sqrt(k2);
	Float mu2= (kz*kz)/k2;
	
	size_t index= nckz*(nc*ix + iy) + kz;
	Float d_re= delta_k[index][0];
	Float d_im= delta_k[index][1];

	// Legendre polynomial
	// P_l = (2 l + 1)/2*int_0^1 P(k) P_l(mu) dmu
	// P_2 = (3 mu^2 + 1)/2
	// P_4 = (35 mu^4 - 30 mu^2 + 3)/8
	Float l_2= 7.5*mu2 - 2.5;
	Float l_4= 1.125*(35.0*mu2*mu2 - 30.0*mu2 + 3.0);
	Float delta2= pk_fac*(d_re*d_re + d_im*d_im);

	if(k < iknq) {
	ps->add(k,
		delta2*corr_xyz - shot_noise,
		l_2*(delta2*corr_xyz - shot_noise),
		l_4*(delta2*corr_xyz - shot_noise));
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
	float l_2= 7.5*mu2 - 2.5;
	float l_4= 1.125*(35.0f*mu2*mu2 - 30.0f*mu2 + 3.0);

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
