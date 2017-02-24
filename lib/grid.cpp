#include <chrono>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "msg.h"
#include "error.h"
#include "grid.h"

static double time_init= 0.0;
static double time_fft= 0.0;

void grid_print_time()
{
  msg_printf(msg_verbose, "Time grid_init %le\n", time_init);
  msg_printf(msg_verbose, "Time grid_fft %le\n",  time_fft);
}

//
// class Grid
//
Grid::Grid(const int nc_) :
  nc(nc_), ncz(2*(nc_/2 + 1)), mode(grid_real_space), boxsize(0.0),
  total_weight(0.0), n_mas(0)
{
  auto ts = std::chrono::high_resolution_clock::now();
  
  //const size_t ncz= 2*(nc/2 + 1);
  x0[0] = x0[1]= x0[2]= 0.0;
  
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

  auto te = std::chrono::high_resolution_clock::now();
  time_init += std::chrono::duration<double>(te - ts).count();
}

void Grid::fft()
{
  auto ts = std::chrono::high_resolution_clock::now();
  
  if(mode != grid_real_space) {
    msg_printf(msg_error,
	       "Error: trying to FFT a grid not in real-space mode\n");
    throw AssertionError();
  }
  
  FFTW(execute)(plan);

  mode = grid_fourier_space;
  auto te = std::chrono::high_resolution_clock::now();
  time_fft += std::chrono::duration<double>(te - ts).count();
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
// Grid algorithms
//
void grid_compute_fluctuation(Grid& grid1, const Grid& grid_rand)
{
  // Compute the fluctuation field: grid1 -= grid_rand
  //
  // Input:
  //     grid1 is the grid of data
  //     grid_rand is the grid of randoms
  //
  // Output:
  //     grid1 become a fluctuation field
  //
  
  auto ts = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Compute fluctuation");

  const size_t nc= grid1.nc;
  const size_t ncz= grid1.ncz;
  assert(grid_rand.nc == nc);

  const double alpha = grid1.total_weight / grid_rand.total_weight;

  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      size_t index= (ix*nc + iy)*ncz;
      for(size_t iz=0; iz<nc; ++iz) {
	grid1[index] -= alpha*grid_rand[index];
	index++;
      }
    }
  }

  // time duration
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time fluctuation %le\n",
	     std::chrono::duration<double>(te - ts).count());
}

void grid_compute_fluctuation_homogeneous(Grid& grid1)
{
  // Compute the fluctuation field: data_grid - nbar,
  // when the nbar is homogeneous, i.e., no randoms catalogue necessary
  //
  // Input:
  //     grid1 is the grid of data
  //
  // Output:
  //     grid1 become a fluctuation field
  //
  
  auto ts = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Compute fluctuation without randoms");

  const size_t nc= grid1.nc;
  const size_t ncz= grid1.ncz;
  const double nbar = grid1.total_weight / (nc*nc*nc);

  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      size_t index= (ix*nc + iy)*ncz;
      for(size_t iz=0; iz<nc; ++iz) {
	grid1[index] -= nbar;
	index++;
      }
    }
  }

  // time duration
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time fluctuation homogeneous %le\n",
	     std::chrono::duration<double>(te - ts).count());
}

