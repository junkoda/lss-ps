//
// Brute-force RR pair count for window function multipoles
//
#include "rr.h"


RRMultipoles::RRMultipoles(const Float r_max_, const Float dr_) :
  dr(dr_), r_max(r_max_)
{
  n= round(r_max/dr) + 1;
  rr0.resize(n, 0.0);
  rr1.resize(n, 0.0);
  rr2.resize(n, 0.0);
  rr3.resize(n, 0.0);
  rr4.resize(n, 0.0);
}



