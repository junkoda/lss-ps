#include <cstdlib>
#include <cstring>
#include <cassert>

#include "msg.h"
#include "error.h"
#include "grid.h"

//
// Grid
//
Grid::Grid(const int nc_) :
  nc(nc_), mode(grid_real_space), boxsize(0.0)
{
  const size_t ncz= 2*(nc/2 + 1);
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
  
  msg_printf(msg_info,
	    "Allocated a grid nc= %d; using %.1lf Mbytes\n",
	       nc, sizeof(Float)*ngrid/(1000.0*1000.0));

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
  const size_t ncz= 2*(nc/2 + 1);
  memset(fx, 0, sizeof(Float)*(nc*nc*ncz));
}

void Grid::write(const char filename[])
{
  FILE* fp= fopen(filename, "w");
  if(fp == 0) {
    msg_printf(msg_error,
	       "Error: unable to write grid data to %s\n", filename);
    throw IOError();
  }

  const int inc= (int) nc;
  const int ncz= 2*(nc/2 + 1);
  
  fwrite(&inc, sizeof(int), 1, fp);
  for(size_t ix= 0; ix<nc; ++ix) {
   for(size_t iy= 0; iy<nc; ++iy) {
    for(size_t iz= 0; iz<nc; ++iz) {
      float f= (float) fx[(ix*nc + iy)*ncz + iz];
      fwrite(&f, sizeof(float), 1, fp);
    }
   }
  }

  fwrite(&inc, sizeof(int), 1, fp);
  fclose(fp);

  msg_printf(msg_info, "Grid data written: %s\n", filename);
}

//
// GridComplex (complex in both real and Fourier space)
//
GridComplex::GridComplex(const int nc_) :
  nc(nc_), mode(grid_real_space), boxsize(0.0)
{
  const size_t ngrid= nc*nc*nc;
  fx= (Complex*) FFTW(malloc)(sizeof(Complex)*ngrid);
  if(fx == 0) {
    msg_printf(msg_fatal,
    "Error: unable to allocate memory for a complex grid nc= %d; %lu Mbytes\n",
	       nc, sizeof(Complex)*ngrid/(1000*1000));
    throw MemoryError();
  }
  
  fk= fx;
  plan= FFTW(plan_dft_3d)(nc, nc, nc, fx, fk, FFTW_FORWARD, FFTW_ESTIMATE);
  
  msg_printf(msg_info,
	    "Allocated a complex grid nc= %d; using %.1lf Mbytes\n",
	       nc, sizeof(Complex)*ngrid/(1000.0*1000.0));
}

GridComplex::~GridComplex()
{
  FFTW(free)(fx);
  FFTW(destroy_plan)(plan);
}

void GridComplex::clear()
{
  memset(fx, 0, sizeof(Complex)*(nc*nc*nc));
}

void GridComplex::fft()
{
  if(mode != grid_real_space) {
    msg_printf(msg_error,
	       "Error: trying to FFT a grid not in real-space mode\n");
    throw AssertionError();
  }
  
  FFTW(execute)(plan);

  mode = grid_fourier_space;
}
