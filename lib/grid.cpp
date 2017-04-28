#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "msg.h"
#include "error.h"
#include "grid.h"

static double time_init= 0.0;
static double time_fft= 0.0;

using namespace std;

void grid_print_time()
{
  msg_printf(msg_verbose, "Time grid_init %le\n", time_init);
  msg_printf(msg_verbose, "Time grid_fft %le\n",  time_fft);
}

//
// class Grid
//
Grid::Grid(const int nc_) :
  nc(nc_), ncz(2*(nc_/2 + 1)), mode(GridMode::real_space), offset(0), boxsize(0),
  total_weight(0), n_mas(0)
{
  auto ts = std::chrono::high_resolution_clock::now();

  x0_box[0] = x0_box[1]= x0_box[2]= 0;
  
  const size_t ngrid= nc*nc*ncz;
  fx= (Float*) FFTW(malloc)(sizeof(Float)*ngrid);
  if(fx == 0) {
    msg_printf(msg_fatal,
	    "Error: unable to allocate memory for a grid nc= %d; %lu Mbytes\n",
	       nc, sizeof(Float)*ngrid/(1000*1000));
    throw MemoryError();
  }
  
  fk= (complex<Float>*) fx;
  plan= FFTW(plan_dft_r2c_3d)(nc, nc, nc, fx, (FFTW(complex)*) fx,
			      FFTW_ESTIMATE);
  
  msg_printf(msg_info,
	    "Allocated a grid nc= %d; using %.1lf Mbytes\n",
	       nc, sizeof(Float)*ngrid/(1000.0*1000.0));

  auto te = std::chrono::high_resolution_clock::now();
  time_init += std::chrono::duration<double>(te - ts).count();
}

void Grid::fft()
{
  auto ts = std::chrono::high_resolution_clock::now();

  if(mode != GridMode::real_space) {
    msg_printf(msg_error,
	       "Error: trying to FFT a grid not in real-space mode\n");
    throw AssertionError();
  }
  
  FFTW(execute)(plan);

  mode = GridMode::fourier_space;
  auto te = std::chrono::high_resolution_clock::now();
  time_fft += std::chrono::duration<double>(te - ts).count();
}

Grid::~Grid()
{
  assert(fx);
  FFTW(free)(fx);
  FFTW(destroy_plan)(plan);
}

void Grid::clear()
{
  total_weight= w2_sum= nw2_sum= 0.0;
  np= 0;

  
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
void grid_compute_fluctuation(Grid& grid_data, Grid& grid_rand)
{
  // Compute the fluctuation field: grid_data -= grid_rand
  //
  // Input:
  //     grid_data is the number density grid of data
  //     grid_rand is the number density grid of randoms
  //
  // Output:
  //     grid_data become a fluctuation field
  //
  
  auto ts = std::chrono::high_resolution_clock::now();


  const size_t nc= grid_data.nc;
  const size_t ncz= grid_data.ncz;
  assert(grid_rand.nc == nc);

  const double alpha = grid_data.total_weight / grid_rand.total_weight;
  //cerr << "data total_weight " << grid_data.total_weight << endl;
  //cerr << "rand total_weight " << grid_rand.total_weight << endl;
  
  msg_printf(msg_verbose, "Compute fluctuation, alpha= %e\n", alpha);

  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      size_t index= (ix*nc + iy)*ncz;
      for(size_t iz=0; iz<nc; ++iz) {
	grid_data[index] -= alpha*grid_rand[index];
	index++;
      }
    }
  }

  grid_data.shot_noise = (1.0 + alpha)*grid_rand.w2_sum / grid_rand.nw2_sum;
  grid_rand.shot_noise = 0.0;

  grid_data.pk_normalisation= 1.0/(grid_rand.nw2_sum*alpha);
  //cerr << "normalization= " << grid_rand.nw2_sum*alpha << endl;
  
  grid_rand.pk_normalisation= 1.0/(grid_rand.nw2_sum/alpha);


  // time duration
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time fluctuation %le\n",
	     std::chrono::duration<double>(te - ts).count());
}

void grid_compute_fluctuation_homogeneous(Grid& grid_data)
{
  // Compute the fluctuation field: data_grid - nbar,
  // when the nbar is homogeneous, i.e., no randoms catalogue necessary
  //
  // Input:
  //     grid_data is the number density grid of data
  //
  // Output:
  //     grid_data become a fluctuation field
  //
  
  auto ts = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Compute fluctuation without randoms.\n");

  const size_t nc= grid_data.nc;
  const size_t ncz= grid_data.ncz;
  const double nbar = grid_data.total_weight / (nc*nc*nc);

  for(size_t ix=0; ix<nc; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      size_t index= (ix*nc + iy)*ncz;
      for(size_t iz=0; iz<nc; ++iz) {
	grid_data[index] -= nbar;
	index++;
      }
    }
  }

  const double boxsize= grid_data.boxsize;
  const double vol= boxsize*boxsize*boxsize;
  grid_data.pk_normalisation= vol
                              / (grid_data.total_weight*grid_data.total_weight);
  grid_data.shot_noise = vol / (double) grid_data.np;
  
  // time duration
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time fluctuation homogeneous %le\n",
	     std::chrono::duration<double>(te - ts).count());
}
