#ifndef _CATALOGUE_H
#define _CATALOGUE_H 1

#include "config.h"
#include <vector>

struct Particle {
  Float x[3], v[3];
};

class Catalogue : public std::vector<Particle> {
  
};

//typedef std::vector<Particle> Catalogue;

#endif
