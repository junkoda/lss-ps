#include <cassert>
#include <cmath>
#include <vector>

#include "config.h"
#include "msg.h"
#include "grid.h"
#include "error.h"
#include "py_assert.h"
//#include <iostream>
#include "py_etc.h"

using namespace std;


PyObject* py_etc_gaussian_smoothing(PyObject* self, PyObject* args)
{
  // Apply gaussian smoothing exp(-0.5*(k*sigma)^2) to a grid of white noise
  
  PyObject *py_grid;
  double sigma;
  if(!PyArg_ParseTuple(args, "Od", &py_grid, &sigma)) {
    return NULL;
  }

  Grid* const grid=
    (Grid*) PyCapsule_GetPointer(py_grid, "_Grid");
  py_assert_ptr(grid);

  const int nc= grid->nc;
  const size_t nckz= nc/2 + 1;
  //Float const * const fx= grid->fx;
  const double fac= 2.0*M_PI/grid->boxsize;
  std::complex<Float>* const fk= grid->fk;

  
  vector<double> c1(nc + 1);
  vector<double> gauss(nc + 1);

  for(int i=0; i<nc; ++i) {
    const int ik= i < nc/2 ? i : i - nc;
    const double k= fac*ik;
    const double sinx = sin(M_PI*ik/nc);

    // TSC mass assignment + shot noise correction
    c1[i] = sqrt(1.0 - sinx*sinx + 2.0/15.0*(sinx*sinx*sinx*sinx));

    // Gaussian smoothing
    gauss[i] = exp(-0.5*k*k*sigma*sigma);
  }

  //c1[0] = c1[nc/2] = 1.0;

  //pi*(2*pi/boxsize)*ik/(2*pi/boxsize*nc) = pi*ik/nc

    
  assert(grid->n_mas == 3);
  assert(grid->mode == GridMode::fourier_space);

  for(int ix=0; ix<nc; ++ix) {
    //const int ikx= ix < nc/2 ? ix : ix - nc;
    //const double kx= fac*ikx;
    //const double sinx = sin(M_PI*ikx/nc);
    //const double sn_x = 1.0 - sinx*sinx + 2.0/15.0*(sinx*sinx*sinx*sinx);
    for(int iy=0; iy<nc; ++iy) {
      //const int iky= iy < nc/2 ? iy : iy - nc;
      //const double ky= fac*iky;
      //const double siny = sin(M_PI*iky/nc);
      //const double sn_xy = sn_x*(1.0 - siny*siny + 2.0/15.0*(siny*siny*siny*siny));

      for(int iz=0; iz<nc/2 + 1; ++iz) {
	//const int ikz= iz;
	//const double kz= fac*iz;
	//const double sinz = sin(M_PI*ikz/nc);
	//const double sn_xyz = sn_xy*(1.0 - sinz*sinz + 2.0/15.0*(sinz*sinz*sinz*sinz));

	const double c1_fac= c1[ix]*c1[iy]*c1[iz];
	const double gauss_fac= gauss[ix]*gauss[iy]*gauss[iz];
	
	size_t index= (nc*ix + iy)*nckz + iz;
	fk[index] *= gauss_fac/c1_fac;
	// correct for TSC for shot noise and then apply Gaussian smoothing
      }
    }
  }

  Py_RETURN_NONE;
}
