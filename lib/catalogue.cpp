//
// Catalogue is a collection of particles (galaxy or random particles)
//

#include <cstdio>
#include <vector>
#include "msg.h"
#include "config.h"
#include "catalogue.h"

using namespace std;

void read_ascii_catalogue(const char filename[], Catalogue& cat)
{
  FILE* fp= fopen(filename, "r");
  if(fp != 0) {
    msg_abort("Error: unable to open catalogue file %s\n", filename);
  }

  
  const size_t nbuf= 1024;
  char buf[nbuf];

  Particle p;
  
  while(fgets(buf, nbuf - 1, fp)) {
    if(buf[0] == '#')
      continue;

    
#ifdef DOUBLEPRECISION
    int ret= sscanf(buf, "%le %le %le", p.x, p.x + 1, p.x + 2);
#else
    int ret= sscanf(buf, "%e %e %e", p.x, p.x + 1, p.x + 2);
#endif
    
    if(ret != 3) {
      msg_abort("Error: unable to read 3 numbers, %s\n", buf);
    }

    cat.push_back(p);
  }

  fclose(fp);
}
