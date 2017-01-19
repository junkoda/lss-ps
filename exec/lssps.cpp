//
// Using Boost program options
//   style:   ./options [options] <required arg>
//   example: ./options --x=3 filename
//

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "lssps.h"

using namespace std;
using namespace boost::program_options;

int main(int argc, char* argv[])
{
  //
  // command-line options (Boost program_options)
  //
  options_description opt("options [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "filename")
    ("nc", value<int>()->default_value(64), "integer n")
    ("boxsize", value<double>()->default_value(600.0, "600"), "boxsize")
    ("alias-correction", value<int>()->default_value(1), "0: no correction, 1: Jing et al, 2: interlacing\n")
    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt; 
    return 0;
  }

  const int nc= vm["nc"].as<int>(); assert(nc > 0);

  const double boxsize= vm["boxsize"].as<double>();

  const string filename= vm["filename"].as<string>();

  const int alias_correction= vm["alias-correction"].as<int>();

  msg_set_prefix("# ");
  
  Catalogue* const cat = new Catalogue();
  catalogue_read_text(cat, filename.c_str());

  PowerSpectrum* const ps= new PowerSpectrum(0.0, 1.0, 0.01);

  if(alias_correction == 0) {
    Grid* const grid = new Grid(nc);
    double x0[] = {0.0, 0.0, 0.0};
  
    mass_assignment_cic(cat, x0, boxsize, grid);
    grid->fft();

    power_spectrum_compute_multipoles_raw(grid, true, -1.6, ps);
  }
  else if(alias_correction == 1) {
    Grid* const grid = new Grid(nc);
    double x0[] = {0.0, 0.0, 0.0};
  
    mass_assignment_cic(cat, x0, boxsize, grid);
    grid->fft();

    power_spectrum_compute_multipoles(grid, true, -1.6, ps);
  }
  else if(alias_correction == 2) {
    // Aliasing correction with 'interlacing'
    GridComplex* const grid = new GridComplex(nc);
    double x0[] = {0.0, 0.0, 0.0};
  
    mass_assignment_interlacing_cic(cat, x0, boxsize, grid);
    grid->fft();
    interlacing(grid);
  
    power_spectrum_compute_multipoles_interlacing(grid, true, -1.6, ps);
  }
  else {
    msg_printf(msg_fatal, "Error: unknown alias correction scheme: %d\n",
	       alias_correction);
  }

  for(int i=0; i<ps->n; ++i) {
    if(ps->nmode_hist[i] > 0)
      printf("%e %e\n", ps->k_hist[i], ps->P0_hist[i]);
  }

  //delete ps;
  //delete grid;
  //delete cat;
  
  return 0;
}
