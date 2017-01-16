#include <cstdlib>
#include <cstring>
#include <cassert>

#include "msg.h"
#include "error.h"
#include "grid.h"

Grid::Grid(const int nc_) :
  nc(nc_), ncz(2*(nc_/2+1)), mode(grid_real_space), boxsize(0.0)
{
  const size_t ngrid= nc*nc*ncz;
  fx= (Float*) FFTW(malloc)(sizeof(Float)*ngrid);
  if(fx == 0) {
    msg_printf(msg_fatal,
	    "Error: unable to allocate memory for a grid nc= %d; %lu Mbytes\n",
	       nc, sizeof(Float)*ngrid/(1000*1000));
    throw MemoryError();
  }
  
  fk= (Complex*) fx;
  plan= FFTW(plan_dft_r2c_3d)(nc, nc, nc, fx, fk, FFTW_ESTIMATE);
  clear();
}

void Grid::fft()
{
  if(mode != grid_real_space) {
    msg_printf(msg_error,
	       "Error: trying to FFT a grid not in real-space mode\n");
    throw AssertionError();
  }
  
  FFTW(execute)(plan);

  mode = grid_fourier_space;
}

Grid::~Grid()
{
  FFTW(free)(fx);
  FFTW(destroy_plan)(plan);
}

void Grid::clear()
{
  memset(fx, 0, sizeof(Float)*(nc*nc*ncz));
}
