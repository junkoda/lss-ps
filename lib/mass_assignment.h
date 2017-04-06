#ifndef MASS_ASSIGNMENT_H
#define MASS_ASSIGNMENT_H 1

#include "catalogue.h"
#include "grid.h"

void mass_assignment(const std::vector<Particle>& cat,
		     const Float x0[], const Float boxsize,
		     const int mas,
		     const bool parallel,
		     Grid* const grid);

#endif
