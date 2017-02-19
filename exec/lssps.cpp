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
    ("shotnoise-correction", value<bool>()->default_value(true), "subtract shot noise true=1/false=0")
    ("alias-correction", value<int>()->default_value(1), "0: no correction, 1: Jing et al, 2: interlacing\n")
    ("write-data-grid", value<string>(), "filename of data grid")
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

  const bool subtract_shotnoise= vm["shotnoise-correction"].as<bool>();
  const int alias_correction= vm["alias-correction"].as<int>();

  msg_set_prefix("# ");
  
  Catalogue* const cat = new Catalogue();
  catalogue_read_text(cat, filename.c_str());

  PowerSpectrum* const ps= new PowerSpectrum(0.0, 1.0, 0.01);

  double x0[] = {0.0, 0.0, 0.0};

  //cerr << "alias correction " << alias_correction << endl;
  
  if(alias_correction == 0) {
    // No correction to power spectrum
    Grid* const grid = new Grid(nc);
  
    mass_assignment_cic(cat, x0, boxsize, grid);

    if(vm.count("write-data-grid")) {
      string ofilename= vm["write-data-grid"].as<string>();
      grid->write(ofilename.c_str());
    }
    
    grid->fft();

    power_spectrum_compute_multipoles_raw(grid, ps);
  }
  else if(alias_correction == 1) {
    Grid* const grid = new Grid(nc);
    
    mass_assignment_cic(cat, x0, boxsize, grid);
    grid->fft();

    power_spectrum_compute_multipoles(grid, ps, subtract_shotnoise, -1.6);
  }
  else if(alias_correction == 2) {
    // Aliasing correction with 'interlacing' using two real grids

    Grid* const grid1 = new Grid(nc);
    mass_assignment_cic(cat, x0, boxsize, grid1);

    Grid* const grid2 = new Grid(nc);
    const double h = 0.5*boxsize/nc;
    //cerr << "h= " << h << endl;

    const double x0_shifted[] = {x0[0] - h, x0[1] - h, x0[2] - h};
    mass_assignment_cic(cat, x0_shifted, boxsize, grid2);


    grid1->fft();
    grid2->fft();
    
    interlacing2(grid1, grid2);
    
    power_spectrum_compute_multipoles_interlacing2(grid1, ps, subtract_shotnoise);
  }
  else if(alias_correction == 3) {
    // Alias correction with 'interlacing' using one complex grid
    GridComplex* const grid = new GridComplex(nc);
    mass_assignment_interlacing_cic(cat, x0, boxsize, grid);
    grid->fft();
    power_spectrum_compute_multipoles_interlacing(grid, ps, subtract_shotnoise);
  }
  else {
    msg_printf(msg_fatal, "Error: unknown alias correction scheme: %d\n",
	       alias_correction);
  }

  for(int i=0; i<ps->n; ++i) {
    if(ps->nmode_hist[i] > 0)
      printf("%e %e %d\n",
	     ps->k_hist[i], ps->P0_hist[i], ps->nmode_hist[i]);
  }

  //delete ps;
  //delete grid;
  //delete cat;
  
  return 0;
}
