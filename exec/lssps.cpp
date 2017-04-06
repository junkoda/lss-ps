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
    ("data", value<string>(), "data filename")
    ("rand", value<string>(), "rand filename")
    ("nc", value<int>()->default_value(64), "integer n")
    ("boxsize", value<double>(), "boxsize")
    ("centre", value<vector<double>>()->multitoken(), "centre of the catalogue")
    ("shotnoise-correction", value<bool>()->default_value(true), "subtract shot noise true=1/false=0")
    ("alias-correction", value<int>()->default_value(1), "0: no correction, 1: Jing et al, 2: interlacing\n")
    ("mas-correction", value<bool>()->default_value(1), "mass assignment window correction\n")
    ("write-data-grid", value<string>(), "filename of data grid")
    ("xyz", value<vector<int>>()->multitoken(), "xyz columns in input file e.g. --xyz=1,2,3")
    ("radec", value<vector<int>>()->multitoken(), "Ra-Dec columns in input file e.g. --radec=1,2")
    ("r", value<int>(), "radius r column in input file e.g. --r=3")
    ("nbar", value<int>()->default_value(0), "mean density column in input file, 0 if no column for nbar")
    ("weights", value<vector<int>>()->multitoken(), "weight columns in input file")
    ("Pest", value<double>()->default_value(0), "set > 0 to use FKP weighting")
    ;
  
  positional_options_description p;
  
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

  // Boxsize
  if(vm.count("rand") == false && vm.count("boxsize") == false) {
    cerr << "Error: --boxsize is necessary when no --rand is given\n";
    return 1;
  }
  double boxsize= vm.count("boxsize") ? vm["boxsize"].as<double>() : 0.0;

  // centre
  vector<double> centre;
  if(vm.count("centre")) {
    centre = vm["centre"].as<vector<double>>();
    if(centre.size() != 3) {
      cerr << "Error: 3 numbers are required for --centre=x y z\n";
      return 1;
    }
  }

  // FKP weight
  const double Pest= vm["Pest"].as<double>();
  const int nc= vm["nc"].as<int>(); assert(nc > 0);

  msg_set_prefix("# ");

  //
  // Read data and randoms from files
  //
  const string data_filename= vm["data"].as<string>();
  Catalogue* const data = new Catalogue();

  string rand_filename;
  const bool use_rand = vm.count("rand");
  Catalogue* const rand = use_rand ? new Catalogue() : 0;
  Grid* const grid_rand = use_rand ? new Grid(nc) : 0;
  if(rand)
    rand_filename= vm["rand"].as<string>();

  
  //catalogue_read_text(data, data_filename.c_str());
  const int inbar= vm["nbar"].as<int>();
  vector<int> iweights;
  if(vm.count("weights"))
    iweights= vm["weights"].as<vector<int>>();
  
  if(vm.count("xyz")) {
    if(vm.count("radec")) {
      msg_printf(msg_error, "Error: cannot specify --radec with --xyz\n");
      throw RuntimeError();
    }
    if(vm.count("r")) {
      msg_printf(msg_error, "Error: cannot specify radius --r with --xyz\n");
      throw RuntimeError();
    }
    vector<int> ixyz= vm["xyz"].as<vector<int>>();
    if(ixyz.size() != 3) {
      msg_printf(msg_error, "Error: three numbers are required for --xyz, e.g., --xyz=1,2,3\n");
      throw RuntimeError();
    }

    catalogue_read(data, data_filename.c_str(), ixyz, iweights, inbar,
		   Pest, XYZ());
    if(rand)
      catalogue_read(rand, rand_filename.c_str(), ixyz, iweights, inbar,
		     Pest, XYZ());
  }
  else if(vm.count("radec") && vm.count("r")) {
    vector<int> ipos= vm["radec"].as<vector<int>>();
    if(ipos.size() != 2) {
      msg_printf(msg_error, "Error: 2 numbers are required for --radec, e.g., --radec=1,2\n");
      throw RuntimeError();
    }

    const int ir= vm["r"].as<int>();
    ipos.push_back(ir);

    catalogue_read(data, data_filename.c_str(), ipos, iweights, inbar,
		   Pest, Spherical(AngleUnit::deg));

    if(rand)
      catalogue_read(rand, rand_filename.c_str(), ipos, iweights, inbar,
		     Pest, Spherical(AngleUnit::deg));
  }
  else {
    msg_printf(msg_error, "Error: not enough information is given to read file; set --xyz=1,2,3 or --radec=1,2 --r=3\n");
  }


  //
  // Compute the dimension of the catalogue if necessary
  //

  // If there is no --random, the corner is (0,0,0) and --boxsize must be given
  double x0[]= {0.0, 0.0, 0.0};

  if(rand) {
    if(vm.count("boxsize") == false || vm.count("centre") == false) {
      double boxsize_computed;
      catalogue_compute_range(*rand, x0, boxsize_computed);

      if(boxsize == 0.0)
	boxsize= boxsize_computed;
    }
  }
  
  if(!centre.empty()) {
    for(int k=0; k<3; ++k)
      x0[k]= centre[k] - 0.5*boxsize;
  }

  msg_printf(msg_info, "Catalogue x0 = %e %e %e, boxsize = %e\n",
	     x0[0], x0[1], x0[2], boxsize);

  //
  // Mass assignment
  //
  
  Grid* const grid_data = new Grid(nc); assert(grid_data);
  //mass_assignment_cic(data, x0, boxsize, false, 0.0, grid1);
  mass_assignment(data, CIC(), x0, boxsize, *grid_data);

  // Write grid data if requested
  if(vm.count("write-data-grid")) {
    string ofilename= vm["write-data-grid"].as<string>();
    grid_data->write(ofilename.c_str());
  }

  if(rand) {
    //catalogue_read_text(rand, rand_filename.c_str());
    assert(grid_rand);
    mass_assignment(rand, CIC(), x0, boxsize, *grid_rand);

    grid_compute_fluctuation(*grid_data, *grid_rand);
  }
  else {
    grid_compute_fluctuation_homogeneous(*grid_data);
  }

  //
  // Shifted grid for interlacing
  //
  const int alias_correction= vm["alias-correction"].as<int>();
  Grid* const grid_data_shifted = (alias_correction == 2) ? new Grid(nc) : 0;
  Grid* const grid_rand_shifted = (alias_correction == 2 && rand) ?
                                  new Grid(nc) : 0;

  if(grid_data_shifted) {
    // TODO This part should be more compact
    const double h = 0.5*boxsize/nc;
    const double x0_shifted[] = {x0[0] - h, x0[1] - h, x0[2] - h};

    mass_assignment(data, CIC(), x0_shifted, boxsize, *grid_data_shifted);

    if(rand) {
      mass_assignment(rand, CIC(), x0_shifted, boxsize, *grid_rand_shifted);
      grid_compute_fluctuation(*grid_data_shifted, *grid_rand_shifted);
    }
    else {
      grid_compute_fluctuation_homogeneous(*grid_data_shifted);
    }

    //interlacing(grid_data_shifted, grid_rand_shifted);
  }

  //
  // FFT real -> Fourier space
  //
  grid_data->fft();
  if(grid_data_shifted)
    grid_data_shifted->fft();

  grid_print_time();
  
  if(grid_data_shifted)
    interlacing(grid_data, grid_data_shifted);
  
  //
  // Compute power spectrum
  //
  const bool subtract_shotnoise= vm["shotnoise-correction"].as<bool>();
  const bool mas_correction = vm["mas-correction"].as<bool>();

  PowerSpectrum* const ps= new PowerSpectrum(0.0, 1.0, 0.01, boxsize);

  power_spectrum_compute_multipoles(grid_data, ps, subtract_shotnoise,
				    mas_correction);
  
  //power_spectrum_compute_multipoles(grid1, ps, subtract_shotnoise, -1.6);


    //msg_printf(msg_fatal, "Error: unknown alias correction scheme: %d\n",
    //alias_correction);

  for(int i=0; i<ps->n; ++i) {
    if(ps->n_modes(i) > 0)
      printf("%e %e %d\n",
	     ps->k(i), ps->P0(i), ps->n_modes(i));
  }

  //delete ps;
  //delete grid;
  //delete cat;

  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time total %le\n",
	     std::chrono::duration<double>(te - ts).count());

  return 0;
}
