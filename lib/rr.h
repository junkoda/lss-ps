#ifndef RR_H
#define RR_H 1

#include <vector>
#include "config.h"

struct RRMultipoles {
  RRMultipoles(const Float r_max_, const Float dr_);
  int n;
  Float dr, r_max;
  std::vector<double> rr0, rr1, rr2, rr3, rr4;
};

//
// Temp
//
template<typename float_type>
void rr_multipoles(float_type const * x,
		   const size_t x_stride,
		   float_type const * weight,
		   const size_t weight_stride,
		   float_type const * nbar,
		   const size_t nbar_stride,
		   const size_t np,
		   RRMultipoles* const rr)
{
  double w= 1.0;
  const int nbin= rr->n;
  const Float dr= rr->dr;
  long double normalisation= 0;
  
  for(size_t i=0; i<np; ++i) {
    float_type w_i = weight == nullptr ? 1.0 : *weight;
    float_type nb = nbar == nullptr ? 1.0 : *nbar;
    float_type const * weight_j= weight;
    float_type const * x_j= x;

    normalisation += nb*w_i*w_i;

    for(size_t j=i+1; j<np; ++j) {
      x_j = (float_type*) ((char*) x_j + x_stride);
      if(weight_j) {
	weight_j = (float_type*) ((char*) weight_j + weight_stride);
	w = w_i*(*weight_j);
      }
      
      Float r[3];
      r[0]= x_j[0] - x[0];
      r[1]= x_j[1] - x[1];
      r[2]= x_j[2] - x[2];

      Float x_norm= sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
      Float r_norm= sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
      
      Float nu= (x[0]*r[0] + x[1]*r[1] + x[2]*r[2])/(x_norm*r_norm);

      Float nu2= nu*nu;
      Float P2= 1.5*nu2 - 0.5;
      Float P3= 2.5*nu2*nu - 1.5*nu;
      Float P4= 0.125*(35.0*nu2*nu2 - 30.0*nu2 + 3.0);

      // exchange i <-> j
      Float x_norm_j= sqrt(x_j[0]*x_j[0] + x_j[1]*x_j[1] + x_j[2]*x_j[2]);
      Float nu_j= -(x_j[0]*r[0] + x_j[1]*r[1] + x_j[2]*r[2])/(x_norm_j*r_norm);
      Float nu2_j= nu_j*nu_j;
      Float P2_j= 1.5*nu2_j - 0.5;
      Float P3_j= 2.5*nu2_j*nu_j - 1.5*nu_j;
      Float P4_j= 0.125*(35.0*nu2_j*nu2_j - 30.0*nu2_j + 3.0);


      int ibin= static_cast<int>(r_norm/dr);
      if(ibin < nbin) {
	rr->rr0[ibin] += 2.0*w;
	rr->rr1[ibin] += w*(nu + nu_j);
	rr->rr2[ibin] += w*(P2 + P2_j);
	rr->rr3[ibin] += w*(P3 + P3_j);
	rr->rr4[ibin] += w*(P4 + P4_j);
      }
    }
    
    x = (float_type*) ((char*) x + x_stride);
    if(weight)
      weight = (float_type*) ((char*) weight + weight_stride);
    if(nbar)
      nbar   = (float_type*) ((char*) nbar   + nbar_stride);
  }

  
  for(int i=0; i<nbin; ++i) {
    double r_left= i*dr;
    double r_right= (i + 1.0)*dr;
    double vol= 4.0*M_PI/3.0*(pow(r_right, 3) - pow(r_left, 3));

    double fac= 1.0/normalisation/vol;
    
    rr->rr0[i] = fac*rr->rr0[i];
    rr->rr1[i] = 3.0*fac*rr->rr1[i];
    rr->rr2[i] = 5.0*fac*rr->rr2[i];
    rr->rr3[i] = 7.0*fac*rr->rr3[i];
    rr->rr4[i] = 9.0*fac*rr->rr4[i];
  }
}

#endif
