#ifndef _CATALOGUE_H
#define _CATALOGUE_H 1

#include "config.h"
#include <vector>

struct Particle {
  Float x[3], w;
};

class Catalogue : public std::vector<Particle> {
  
};

void catalogue_read_text(Catalogue* const cat, const char filename[]);
  
#endif
