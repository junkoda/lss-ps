//
// Using Boost program options
//   style:   ./options [options] <required arg>
//   example: ./options --x=3 filename
//

#include <iostream>
#include <string>
#include <chrono>
#include <boost/program_options.hpp>

#include "lssps.h"


using namespace std;
using namespace boost::program_options;

int main(int argc, char* argv[])
{
  auto ts = std::chrono::high_resolution_clock::now();
  
  //
  // command-line options (Boost program_options)
  //
  options_description opt("options [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    //("filename,f", value<string>(), "filename")
    ("data", value<string>(), "data filename")
    ("rand", value<string>(), "rand filename")
    ("nc", value<int>()->default_value(64), "integer n")
    ("boxsize", value<double>()->default_value(600.0, "600"), "boxsize")
    ("shotnoise-correction", value<bool>()->default_value(true), "subtract shot noise true=1/false=0")
    ("alias-correction", value<int>()->default_value(1), "0: no correction, 1: Jing et al, 2: interlacing\n")
    ("write-data-grid", value<string>(), "filename of data grid")
    ;
  
  positional_options_description p;
  //p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help")) {
    cout << opt; 
    return 0;
  }

  if(!vm.count("data")) {
    cerr << "Error: --data file name not given\n";
    return 1;
  }

  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  const double boxsize= vm["boxsize"].as<double>();

  //const string filename= vm["filename"].as<string>();

  msg_set_prefix("# ");

  // Read data and create data grid
  const string data_filename= vm["data"].as<string>();
  Catalogue* const data = new Catalogue();
  const double x0[]= {0.0, 0.0, 0.0};

  catalogue_read_text(data, data_filename.c_str());

  Grid* const grid1 = new Grid(nc); assert(grid1);
  mass_assignment(data, CIC(), x0, boxsize, false, 0.0, *grid1);

  // Write grid data
  if(vm.count("write-data-grid")) {
    string ofilename= vm["write-data-grid"].as<string>();
    grid1->write(ofilename.c_str());
  }

  // Read randoms and convert grid1 to fluctuation, data - random.
  const bool use_rand = vm.count("rand");
  Catalogue* const rand = use_rand ? new Catalogue() : 0;
  Grid* const grid2 = use_rand ? new Grid(nc) : 0;

  if(rand) {
    const string rand_filename= vm["rand"].as<string>();
    catalogue_read_text(rand, rand_filename.c_str());
    assert(grid2);

    mass_assignment(rand, CIC(), x0, boxsize, false, 0.0, *grid2);

    grid_compute_fluctuation(*grid1, *grid2);
  }
  else {
    grid_compute_fluctuation_homogeneous(*grid1);
  }

  //
  // Shifted grid for interlacing
  //
  const int alias_correction= vm["alias-correction"].as<int>();
  Grid* const grid1_shifted = (alias_correction == 2) ? new Grid(nc) : 0;

  if(grid1_shifted) {
    // TODO This part should be more compate
    const double h = 0.5*boxsize/nc;
    const double x0_shifted[] = {x0[0] - h, x0[1] - h, x0[2] - h};

    mass_assignment(data, CIC(), x0_shifted, boxsize, false, 0.0,
		    *grid1_shifted);

    if(rand) {
      mass_assignment(rand, CIC(), x0_shifted, boxsize, false, 0.0, *grid2);
      grid_compute_fluctuation(*grid1_shifted, *grid2);
    }
    else {
      grid_compute_fluctuation_homogeneous(*grid1_shifted);
    }

    interlacing(grid1, grid1_shifted);
  }

  //
  // FFT real -> Fourier space
  //
  grid1->fft();
  if(grid1_shifted)
    grid1_shifted->fft();
  
  //
  // Compute power spectrum
  //
  const bool subtract_shotnoise= vm["shotnoise-correction"].as<bool>();

  PowerSpectrum* const ps= new PowerSpectrum(0.0, 1.0, 0.01);

  power_spectrum_compute_multipoles(grid1, ps, subtract_shotnoise);
  
  //power_spectrum_compute_multipoles(grid1, ps, subtract_shotnoise, -1.6);


    //msg_printf(msg_fatal, "Error: unknown alias correction scheme: %d\n",
    //alias_correction);

  for(int i=0; i<ps->n; ++i) {
    if(ps->nmode_hist[i] > 0)
      printf("%e %e %d\n",
	     ps->k_hist[i], ps->P0_hist[i], ps->nmode_hist[i]);
  }

  //delete ps;
  //delete grid;
  //delete cat;

  grid_print_time();
  
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time total %le\n",
	     std::chrono::duration<double>(te - ts).count());

  return 0;
}
