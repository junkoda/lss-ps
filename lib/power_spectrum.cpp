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
// class Power Spectrum
//
PowerSpectrum::PowerSpectrum()
{  

}

PowerSpectrum::PowerSpectrum(const int n_) :
  n(n_), nmodes(0.0, n),
  k(0.0, n), p0(0.0, n), p2(0.0, n), p4(0.0, n)
{

}
  
PowerSpectrum& PowerSpectrum::operator+=(const PowerSpectrum& ps)
{
  nmodes += ps.nmodes;
  
  k += ps.k;
  p0 += ps.p0;
  p2 += ps.p2;
  p4 += ps.p4;

  return *this;
}
