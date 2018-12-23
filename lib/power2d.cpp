#include <vector>
#include <valarray>
#include <chrono>
#include <cassert>

#include "config.h"
#include "msg.h"
#include "mas_correction.h"
#include "power2d.h"

using namespace std;

//
// Function objects for compute_power2d_template
//
class CrossPower2dRe {
public:
  explicit CrossPower2dRe(Grid const * const grid1_,
			  Grid const * const grid2_) :
    grid(grid1_), grid2(grid2_) { }

  double operator()(const size_t index) const {
    return grid->fk[index].real()*grid2->fk[index].real()
         + grid->fk[index].imag()*grid2->fk[index].imag();
  }
  Grid const * const grid;
  Grid const * const grid2;
};

class CrossPower2dIm {
public:
  explicit CrossPower2dIm(Grid const * const grid1_,
			  Grid const * const grid2_) :
    grid(grid1_), grid2(grid2_) { }

  double operator()(const size_t index) const {
    return grid->fk[index].imag()*grid2->fk[index].real()
         - grid->fk[index].real()*grid2->fk[index].imag();
  }
  Grid const * const grid;
  Grid const * const grid2;
};


//
// Class PowerSpectrum2D
//
class PowerSpectrum2D {
 public:
  PowerSpectrum2D();
  PowerSpectrum2D(const int n);
  PowerSpectrum2D& operator+=(const PowerSpectrum2D& ps);
  
  PowerSpectrum2D(const PowerSpectrum2D&) = delete;
  PowerSpectrum2D& operator=(const PowerSpectrum2D& ps) = delete;

  valarray<double> nmodes;
  valarray<double> k, mu, p2d;
};

PowerSpectrum2D::PowerSpectrum2D()
{

}

PowerSpectrum2D::PowerSpectrum2D(const int n) :
  nmodes(0.0, n), k(0.0, n), mu(0.0, n), p2d(0.0, n)
{

}

PowerSpectrum2D& PowerSpectrum2D::operator+=(const PowerSpectrum2D& ps)
{
  nmodes += ps.nmodes;
  k += ps.k;
  mu += ps.mu;
  p2d += ps.p2d;

  return *this;
}


//
// 2D power computation algorithm
//
template <typename F>
void compute_power2d_template(const double k_min,
			      const double dk,
			      const int nk, const int nmu,
			      const F f,
			      const double shot_noise,
			      const bool correct_mas,
			      const int los,
			      double * const nmodes_out,
			      double * const k_out,
			      double * const mu_out,
			      double * const P_out
			      )
{
  auto ts = std::chrono::high_resolution_clock::now();

  assert(0 <= los && los < 3);
  assert(nmu > 0);

  Grid const * const grid = f.grid; assert(grid);
  const int nc= grid->nc;
  const double boxsize= grid->boxsize;
  const int n_mas= correct_mas ? grid->n_mas : 0;

  const vector<Float>& mas_correction_array=
    mas_correction_get_array(nc, n_mas);

  assert(nc > 0);
  assert(boxsize > 0.0);
  assert(grid->mode == GridMode::fourier_space);
  assert(mas_correction_array.size() == (size_t) nc);

  // binning
  const double k_fundamental= 2.0*M_PI/boxsize;
  //const double k_nq= M_PI/boxsize*nc;
  const double ik_min= k_min/k_fundamental;
  const double ik_max= (k_min + nk*dk)/k_fundamental;
  const double idk= dk/k_fundamental;

  PowerSpectrum2D P(nk*nmu);

  const size_t nckz= nc/2+1;
  const int ik_nq= nc/2;
  const int ik_max2 = ik_max*ik_max;

  msg_printf(msg_info, "Shot noise subtraction: %e\n", shot_noise);
  msg_printf(msg_info, "Pk normalisation %e\n", grid->pk_normalisation);

#ifdef _OPENMP
  #pragma omp parallel num_threads(omp_get_max_threads())
#endif
  {
    PowerSpectrum2D ps_local(nk*nmu);

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
	  int ik= static_cast<int>(floor((kmag - ik_min)/idk));

	  double mu= k[los]/kmag;
	  int imu= static_cast<int>(floor(mu*nmu));
	  if(imu == nmu) imu= nmu - 1;
	
	  if(0 <= ik && ik < nk) {
	    size_t index= nckz*(nc*ix + iy) + iz;
	    int ibin = ik*nmu + imu;
	    double corr_xyz = corr_xy * mas_correction_array[iz];

	    ps_local.nmodes[ibin] += 1.0;
	    ps_local.k[ibin] += kmag;
	    ps_local.mu[ibin] += mu;
	    ps_local.p2d[ibin] += f(index)*corr_xyz;
	  }	
	}
      }
    }
			   

    // Local power spectrum in each core is added up to total
#ifdef _OPENMP
    #pragma omp critical (__COLLECT_PK_2D__)
#endif
    {
      P += ps_local;
    }
  } // end of omp parallel

  // Get the mean by dividing by the number of modes
  const double pk_fac= grid->pk_normalisation; 
  for(int ik=0; ik<nk; ++ik) {
    for(int imu=0; imu<nmu; ++imu) {
      int i= ik*nmu + imu;
      if(P.nmodes[i] > 0) {
	k_out[i]= k_fundamental*P.k[i]/P.nmodes[i];
	mu_out[i]= P.mu[i]/P.nmodes[i]; 
      
	P_out[i] = pk_fac*P.p2d[i]/P.nmodes[i] - shot_noise;
      }
      else {
	k_out[i]= 0.0;
	mu_out[i]= 0.0;
	P_out[i]= 0.0;
      }
    }
  }
  
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_info, "Time power2d %e\n",
             std::chrono::duration<double>(te - ts).count());
}


void power2d_compute(Grid const * const grid1,
		     Grid const * const grid2,
		     const int real_imag,
		     const double k_min,
		     const double dk,
		     const int nk, const int nmu,
		     const double shot_noise,
		     const bool correct_mas,
		     const int los,
		     double * const nmodes_out,
		     double * const k_out,
		     double * const mu_out,
		     double * const P_out)
{
  if(real_imag == 0) {
    compute_power2d_template(k_min, dk, nk, nmu,
			     CrossPower2dRe(grid1, grid2),
			     shot_noise, correct_mas, los,
			     nmodes_out, k_out, mu_out, P_out);
  }
  else if(real_imag == 1) {
    compute_power2d_template(k_min, dk, nk, nmu,
			     CrossPower2dIm(grid1, grid2),
			     shot_noise, correct_mas, los,
			     nmodes_out, k_out, mu_out, P_out);
  }
  else {
    msg_abort("Unknown real_imag value %d; must be 0 or 1", real_imag);
  }
}
