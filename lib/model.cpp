#include <complex>
#include <cstdio>
#include <cassert>
#include <vector>
#include <cmath>

#include "config.h"
#include "msg.h"
#include "error.h"
#include "cosmology.h"
#include "grid.h"
#include "model.h"

using namespace std;

static void read_power_spectrum(const char filename[],
				vector<double>& v_k, vector<double>& v_P,
				const double fac)
{
  char buf[256];
  
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal,
	       "Error: Unable to open power spectrum file for model. %s\n",
	       filename);
    throw IOError();
  }

  while(fgets(buf, 255, fp)) {
    if(buf[0] == '#')
      continue;
    double k, P;
    
    int ret= sscanf(buf, "%le %le", &k, &P);
    if(ret != 2) {
      msg_printf(msg_fatal,
		 "Error: Unable to read k P from power spectrum file: %s\n",
		 filename);
      msg_printf(msg_fatal, "line: %s", buf);
      throw IOError();      
    }

    v_k.push_back(k);
    v_P.push_back(fac*P);
  }

  fclose(fp);

  msg_printf(msg_verbose,
	     "Read %d lines from power spectrum file for a model\n",
	     (int) v_k.size());

  assert(v_k.size() == v_P.size());
}

//
// ModelLinear
//

//
// ModelLinear class
//
ModelLinear::ModelLinear(const char filename[],
			 const double omega_m,
			 const double redshift,
			 const double bias,
			 const bool redshift_space,
			 const double sigma_v_damping) :
  b(bias), sigma_v(sigma_v_damping)
{
  const double a= 1.0/(1.0 + redshift);
  double growth= 1.0;
  if(redshift != 0.0) {
    growth= cosmology_growth_D(a, omega_m);
  }
  const double fac= growth*growth;

  if(redshift_space) {
    f= cosmology_growth_f(a, omega_m);
  }
  else {
    f= 0.0;
  }

  
  vector<double> v_k, v_P;
  read_power_spectrum(filename, v_k, v_P, fac);

  spline= gsl_spline_alloc(gsl_interp_cspline, v_k.size());
  acc= gsl_interp_accel_alloc();
  gsl_spline_init(spline, v_k.data(), v_P.data(), v_k.size());

  k_min = v_k.front();
  k_max = v_k.back();
}

ModelLinear::~ModelLinear()
{
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

double ModelLinear::operator()(const double k, const double mu) const
{
  if(k < k_min || k > k_max)
    return 0.0;

  const double pk= gsl_spline_eval(spline, k, acc);
  const double fac= (b + f*mu*mu)*exp(-k*sigma_v);
  
  return fac*fac*pk;
}

//
// Discrete power spectrum multipoles 
//
vector<Float> mas_window_array(const int nc,
			       const int mas_nc, const int mas_deg)
{
  vector<Float> v(nc, 1.0);

  if(mas_deg == 0)
    return v;
  
  const int knq = nc/2;
  const Float fac= M_PI/mas_nc;

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i=1; i<nc; ++i) {
    int k= i <= knq ? i : i - nc;
    Float sinc = sin(fac*k)/(fac*k);
    v[i] = pow(sinc, mas_deg);
  }
  
  return v;
}

//
// Algorithms
//
PowerSpectrum* model_compute_discrete_multipoles(const int nc,
						 const Float boxsize,
						 const double k_min,
						 const double k_max,
						 const double dk,
						 const Model& f,
						 const int mas_window_nc,
						 const int mas_window_deg)
{
  const std::vector<Float> vmas=
    mas_window_array(nc, mas_window_nc, mas_window_deg);

  const double fac= 2.0*M_PI/boxsize;
  const int n= (int) round((k_max - k_min)/dk);
  PowerSpectrum* const ps= new PowerSpectrum(n);

  for(int ix=0; ix<nc; ++ix) {
    double kx= ix <= nc/2 ? fac*ix : fac*(ix - nc);
    double corr_x= vmas[ix];

    for(int iy=0; iy<nc; ++iy) {
      double ky= iy <= nc/2 ? fac*iy : fac*(iy - nc);
      double corr_xy= corr_x*vmas[iy];

      int iz0 = !(kx > 0.0 || (kx == 0.0 && ky > 0.0));

      // Avoid double counting on kz=plain
      // k=(0,0,0) dropped because this is 0
      // iz0= 0 if kx>0 or (kx == 0 and ky > 0)
      //      1 otherwize
      //
      
      for(int iz=iz0; iz<nc/2+1; ++iz) {
	double kz= fac*iz;
	double corr_xyz= corr_xy*vmas[iz];
      
	double k= sqrt(kx*kx + ky*ky + kz*kz);

	int i= (int) floor((k - k_min)/dk);
      
	if(0 <= i && i < n) {
	  double mu= kz/k;
	  double mu2 = mu*mu;
	  double l2= 7.5*mu2 - 2.5;
	  double l4= (1.125*35.0)*mu2*mu2 - (1.125*30.0)*mu2 + (1.125*3.0);
	  double P= f(k, mu)*corr_xyz;
	  
	  ps->nmodes[i] += 1.0;
	  ps->k[i] += k;	  
	  ps->p0[i] += P;
	  ps->p2[i] += l2*P;
	  ps->p4[i] += l4*P;
	}
      }  
    }
  }

  ps->k  /= ps->nmodes;
  ps->p0 /= ps->nmodes;
  ps->p2 /= ps->nmodes;
  ps->p4 /= ps->nmodes;

  return ps;
};

//
// Apply 3D window function
//
double model_apply_window_3d(const Model& model,
			     Grid const * const grid,
			     const double k)
{
  // Compute convolution
  // P~(k) = int d3k'/(2pi)^3 |W(k')|^2 P(|k-k'|)
  // for kvec = (k, 0, 0)
  //
  // grid is |W(k)|^2
  const int nc= grid->nc;
  const int nckz= nc/2 + 1;
  const double iknq= nc/2;

  const double fac= 2.0*M_PI/grid->boxsize;
  //const double knq= M_PI/grid->boxsize*nc;

  complex<Float> const * const fk= grid->fk;

  long double P_conv= 0;
  
  for(int ix=0; ix<nc; ++ix) {
    const double kx= ix <= iknq ? fac*ix : fac*(ix - nc);
    for(int iy=0; iy<nc; ++iy) {
      const double ky= iy <= iknq ? fac*iy : fac*(iy - nc);

      int iz0 = !(kx > 0.0 || (kx == 0.0 && ky > 0.0));
      for(int iz=iz0; iz<nckz; ++iz) {
	const double kz= fac*iz;
	const double k1x= k - kx;
	const double k1= sqrt(k1x*k1x + ky*ky + kz*kz);

	const size_t index= (ix*static_cast<size_t>(nc) + iy)*nckz + iz;
	const double wk2= fk[index].real(); // |w(k_1)|^2
	P_conv += model(k1, 0.0)*wk2; // P(k + kx, ky, kz)
      }
    }
  }

  // x2 for the other half of the complex d3k'
  P_conv= 2.0*P_conv + model(k, 0.0)*fk[0].real();

  const double vol= pow(grid->boxsize, 3.0);
  // \int d3k/(2pi)^3 = 1/vol sum_k

  return P_conv/vol;
}
