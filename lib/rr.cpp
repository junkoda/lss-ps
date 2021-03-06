//
// Brute-force RR pair count for window function multipoles
//
#include "rr.h"

using namespace std;

RRMultipoles::RRMultipoles(const Float r_max_, const Float dr_) :
  dr(dr_), r_max(r_max_)
{
  n= round(r_max/dr) + 1;

  rr= new vector<double>[(nn + 1)*(nl + 1)];

  for(int i=0; i<=nn; ++i) {
    for(int j=0; j<=nl; ++j) {
      rr[i*nl + j].resize(n, 0.0);
    }
  }
}

RRMultipoles::~RRMultipoles()
{
  delete [] rr;
}




