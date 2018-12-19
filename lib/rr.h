#ifndef RR_H
#define RR_H 1

#include <vector>
#include <cassert>
#include "config.h"

struct RRMultipoles {
  RRMultipoles(const Float r_max_, const Float dr_);
  ~RRMultipoles();
  int n;
  Float dr, r_max;

  std::vector<double>* rr;
  
  // rrl_n -> Q_l^{(n)}
  static const int nn=5, nl=5;
};

//
// RR window function end-point multipoles
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
  // nbar: number density of the *randoms*
  double w= 1.0;
  const int nbin= rr->n;
  const int nl= rr->nl;
  const Float dr= rr->dr;
  long double normalisation= 0;

  assert(rr->nl >= 5);
  assert(rr->nn >= 5);
  for(int i=0; i<rr->nn; ++i) {
    for(int j=0; j<rr->nl; ++j) {
      assert(rr->rr[nl*i + j].size() == static_cast<size_t>(nbin));
    }
  }
  
  // rr[n][l]
  // n=0
  std::vector<double>& rr00= rr->rr[0];
  std::vector<double>& rr01= rr->rr[1];
  std::vector<double>& rr02= rr->rr[2];
  std::vector<double>& rr03= rr->rr[3];
  std::vector<double>& rr04= rr->rr[4];

  // n=1
  std::vector<double>& rr10= rr->rr[1*nl + 0];
  std::vector<double>& rr11= rr->rr[1*nl + 1];
  std::vector<double>& rr12= rr->rr[1*nl + 2];
  std::vector<double>& rr13= rr->rr[1*nl + 3];
  std::vector<double>& rr14= rr->rr[1*nl + 4];

  // n=2
  std::vector<double>& rr20= rr->rr[2*nl + 0];
  std::vector<double>& rr21= rr->rr[2*nl + 1];
  std::vector<double>& rr22= rr->rr[2*nl + 2];
  std::vector<double>& rr23= rr->rr[2*nl + 3];
  std::vector<double>& rr24= rr->rr[2*nl + 4];

  // n=3
  std::vector<double>& rr30= rr->rr[3*nl + 0];
  std::vector<double>& rr31= rr->rr[3*nl + 1];
  std::vector<double>& rr32= rr->rr[3*nl + 2];
  std::vector<double>& rr33= rr->rr[3*nl + 3];
  std::vector<double>& rr34= rr->rr[3*nl + 4];

  // n=4
  std::vector<double>& rr40= rr->rr[4*nl + 0];
  std::vector<double>& rr41= rr->rr[4*nl + 1];
  std::vector<double>& rr42= rr->rr[4*nl + 2];
  std::vector<double>& rr43= rr->rr[4*nl + 3];
  std::vector<double>& rr44= rr->rr[4*nl + 4];

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

      Float r_norm= sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
      
      Float x2= x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
      Float x_norm= sqrt(x2);
      Float x3 = x_norm*x2;
      Float x4 = x2*x2;
      
      Float nu= (x[0]*r[0] + x[1]*r[1] + x[2]*r[2])/(x_norm*r_norm);

      Float nu2= nu*nu;
      Float P2= 1.5*nu2 - 0.5;
      Float P3= 2.5*nu2*nu - 1.5*nu;
      Float P4= 0.125*(35.0*nu2*nu2 - 30.0*nu2 + 3.0);

      // exchange i <-> j
      // y = |x_j|
      Float y2= x_j[0]*x_j[0] + x_j[1]*x_j[1] + x_j[2]*x_j[2];
      Float y= sqrt(y2);
      Float y3= y*y2;
      Float y4= y2*y2;

      
      Float nu_j= -(x_j[0]*r[0] + x_j[1]*r[1] + x_j[2]*r[2])/(y*r_norm);
      Float nu2_j= nu_j*nu_j;
      Float P2_j= 1.5*nu2_j - 0.5;
      Float P3_j= 2.5*nu2_j*nu_j - 1.5*nu_j;
      Float P4_j= 0.125*(35.0*nu2_j*nu2_j - 30.0*nu2_j + 3.0);


      int ibin= static_cast<int>(r_norm/dr);
      if(ibin < nbin) {
	// n = 0
	rr00[ibin] += 2.0*w;
	rr01[ibin] += w*(nu + nu_j);
	rr02[ibin] += w*(P2 + P2_j);
	rr03[ibin] += w*(P3 + P3_j);
	rr04[ibin] += w*(P4 + P4_j);

	// n = 1
	rr10[ibin] += w*(1.0/x_norm + 1.0/y);
	rr11[ibin] += w*(nu/x_norm + nu_j/y);
	rr12[ibin] += w*(P2/x_norm + P2_j/y);
	rr13[ibin] += w*(P3/x_norm + P3_j/y);
	rr14[ibin] += w*(P4/x_norm + P4_j/y);

	// n = 2
	rr20[ibin] += w*(1.0/x2 + 1.0/y2);
	rr21[ibin] += w*(nu/x2 + nu_j/y2);
	rr22[ibin] += w*(P2/x2 + P2_j/y2);
	rr23[ibin] += w*(P3/x2 + P3_j/y2);
	rr24[ibin] += w*(P4/x2 + P4_j/y2);

	// n = 3
	rr30[ibin] += w*(1.0/x3 + 1.0/y3);
	rr31[ibin] += w*(nu/x3 + nu_j/y3);
	rr32[ibin] += w*(P2/x3 + P2_j/y3);
	rr33[ibin] += w*(P3/x3 + P3_j/y3);
	rr34[ibin] += w*(P4/x3 + P4_j/y3);

	// n = 4
	rr40[ibin] += w*(1.0/x4 + 1.0/y4);
	rr41[ibin] += w*(nu/x4 + nu_j/y4);
	rr42[ibin] += w*(P2/x4 + P2_j/y4);
	rr43[ibin] += w*(P3/x4 + P3_j/y4);
	rr44[ibin] += w*(P4/x4 + P4_j/y4);
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

    int nl= rr->nl;
    for(int j=0; j<rr->nn; ++j) {
      rr->rr[nl*j    ][i] *= fac;
      rr->rr[nl*j + 1][i] *= 3.0*fac;
      rr->rr[nl*j + 2][i] *= 5.0*fac;
      rr->rr[nl*j + 3][i] *= 7.0*fac;
      rr->rr[nl*j + 4][i] *= 9.0*fac;
    }
  }
}

//
// RR window function end-point multipoles
//
template<typename float_type>
void rr_multipoles_plane_parallel(float_type const * x,
				  const size_t x_stride,
				  float_type const * weight,
				  const size_t weight_stride,
				  float_type const * nbar,
				  const size_t nbar_stride,
				  const size_t np,
				  const int los,
				  RRMultipoles* const rr)
{
  // nbar: number density of the *randoms*
  double w= 1.0;
  const int nbin= rr->n;
  const int nl= rr->nl;
  const Float dr= rr->dr;
  long double normalisation= 0;

  assert(rr->nl >= 5);
  assert(rr->nn >= 5);
  for(int i=0; i<rr->nn; ++i) {
    for(int j=0; j<rr->nl; ++j) {
      assert(rr->rr[nl*i + j].size() == static_cast<size_t>(nbin));
    }
  }
  
  // rr[n][l]
  // n=0
  std::vector<double>& rr00= rr->rr[0];
  std::vector<double>& rr01= rr->rr[1];
  std::vector<double>& rr02= rr->rr[2];
  std::vector<double>& rr03= rr->rr[3];
  std::vector<double>& rr04= rr->rr[4];

  // n=1
  std::vector<double>& rr10= rr->rr[1*nl + 0];
  std::vector<double>& rr11= rr->rr[1*nl + 1];
  std::vector<double>& rr12= rr->rr[1*nl + 2];
  std::vector<double>& rr13= rr->rr[1*nl + 3];
  std::vector<double>& rr14= rr->rr[1*nl + 4];

  // n=2
  std::vector<double>& rr20= rr->rr[2*nl + 0];
  std::vector<double>& rr21= rr->rr[2*nl + 1];
  std::vector<double>& rr22= rr->rr[2*nl + 2];
  std::vector<double>& rr23= rr->rr[2*nl + 3];
  std::vector<double>& rr24= rr->rr[2*nl + 4];

  // n=3
  std::vector<double>& rr30= rr->rr[3*nl + 0];
  std::vector<double>& rr31= rr->rr[3*nl + 1];
  std::vector<double>& rr32= rr->rr[3*nl + 2];
  std::vector<double>& rr33= rr->rr[3*nl + 3];
  std::vector<double>& rr34= rr->rr[3*nl + 4];

  // n=4
  std::vector<double>& rr40= rr->rr[4*nl + 0];
  std::vector<double>& rr41= rr->rr[4*nl + 1];
  std::vector<double>& rr42= rr->rr[4*nl + 2];
  std::vector<double>& rr43= rr->rr[4*nl + 3];
  std::vector<double>& rr44= rr->rr[4*nl + 4];

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

      Float r_norm= sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
      
      Float x2= x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
      Float x_norm= sqrt(x2);
      Float x3 = x_norm*x2;
      Float x4 = x2*x2;

      Float nu= r[los]/r_norm;

      Float nu2= nu*nu;
      Float P2= 1.5*nu2 - 0.5;
      Float P3= 2.5*nu2*nu - 1.5*nu;
      Float P4= 0.125*(35.0*nu2*nu2 - 30.0*nu2 + 3.0);

      // exchange i <-> j
      // y = |x_j|
      Float y2= x_j[0]*x_j[0] + x_j[1]*x_j[1] + x_j[2]*x_j[2];
      Float y= sqrt(y2);
      Float y3= y*y2;
      Float y4= y2*y2;

      Float nu_j= -nu;
      Float nu2_j= nu_j*nu_j;
      Float P2_j= 1.5*nu2_j - 0.5;
      Float P3_j= 2.5*nu2_j*nu_j - 1.5*nu_j;
      Float P4_j= 0.125*(35.0*nu2_j*nu2_j - 30.0*nu2_j + 3.0);


      int ibin= static_cast<int>(r_norm/dr);
      if(ibin < nbin) {
	// n = 0
	rr00[ibin] += 2.0*w;
	rr01[ibin] += w*(nu + nu_j);
	rr02[ibin] += w*(P2 + P2_j);
	rr03[ibin] += w*(P3 + P3_j);
	rr04[ibin] += w*(P4 + P4_j);

	// n = 1
	rr10[ibin] += w*(1.0/x_norm + 1.0/y);
	rr11[ibin] += w*(nu/x_norm + nu_j/y);
	rr12[ibin] += w*(P2/x_norm + P2_j/y);
	rr13[ibin] += w*(P3/x_norm + P3_j/y);
	rr14[ibin] += w*(P4/x_norm + P4_j/y);

	// n = 2
	rr20[ibin] += w*(1.0/x2 + 1.0/y2);
	rr21[ibin] += w*(nu/x2 + nu_j/y2);
	rr22[ibin] += w*(P2/x2 + P2_j/y2);
	rr23[ibin] += w*(P3/x2 + P3_j/y2);
	rr24[ibin] += w*(P4/x2 + P4_j/y2);

	// n = 3
	rr30[ibin] += w*(1.0/x3 + 1.0/y3);
	rr31[ibin] += w*(nu/x3 + nu_j/y3);
	rr32[ibin] += w*(P2/x3 + P2_j/y3);
	rr33[ibin] += w*(P3/x3 + P3_j/y3);
	rr34[ibin] += w*(P4/x3 + P4_j/y3);

	// n = 4
	rr40[ibin] += w*(1.0/x4 + 1.0/y4);
	rr41[ibin] += w*(nu/x4 + nu_j/y4);
	rr42[ibin] += w*(P2/x4 + P2_j/y4);
	rr43[ibin] += w*(P3/x4 + P3_j/y4);
	rr44[ibin] += w*(P4/x4 + P4_j/y4);
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

    int nl= rr->nl;
    for(int j=0; j<rr->nn; ++j) {
      rr->rr[nl*j    ][i] *= fac;
      rr->rr[nl*j + 1][i] *= 3.0*fac;
      rr->rr[nl*j + 2][i] *= 5.0*fac;
      rr->rr[nl*j + 3][i] *= 7.0*fac;
      rr->rr[nl*j + 4][i] *= 9.0*fac;
    }
  }
}

#endif
