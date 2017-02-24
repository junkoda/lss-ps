//
// Catalogue is a collection of particles (galaxy or random particles)
//

#include <vector>
#include <chrono>
#include <cstdio>

#include "config.h"
#include "msg.h"
#include "error.h"
#include "catalogue.h"

using namespace std;

// Read ascii file to Catalogue
// File format:
//   x y z -- white space separated cartisian corrdinate in 1/h Mpc
//   Lines starting with # is a commnet


void catalogue_read_text(Catalogue* const cat, const char filename[])
{
  auto ts = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Reading catalogue %s\n", filename);
    
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal,
	       "Error: unable to open catalogue file %s\n", filename);
    throw IOError();
  }

  
  const size_t nbuf= 1024;
  char buf[nbuf];

  Particle p;
  p.w = 1;
  
  while(fgets(buf, nbuf - 1, fp)) {
    if(buf[0] == '#')
      continue;

#ifdef DOUBLEPRECISION
    int ret= sscanf(buf, "%le %le %le", p.x, p.x + 1, p.x + 2);
#else
    int ret= sscanf(buf, "%e %e %e", p.x, p.x + 1, p.x + 2);
#endif
    
    if(ret != 3) {
      msg_printf(msg_error, "Error: unable to read 3 numbers, %s\n", buf);
      throw IOError();
    }

    cat->push_back(p);
  }

  int ret = fclose(fp);
  if(ret != 0) {
    msg_printf(msg_error,
	       "Error: error upon closing the catalogue file %s\n", filename);
    throw IOError();
  }

  msg_printf(msg_verbose, "Read %lu particles\n", cat->size());
  auto te = std::chrono::high_resolution_clock::now();
  msg_printf(msg_verbose, "Time read catalogue txt %le\n",
	     std::chrono::duration<double>(te - ts).count());

}
