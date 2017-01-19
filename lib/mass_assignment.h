#ifndef MASS_ASSIGNMENT_H
#define MASS_ASSIGNMENT_H 1

#include "catalogue.h"
#include "grid.h"

void mass_assignment_cic(Catalogue const * const cat,
			 const Float x0[], const Float boxsize,
			 Grid* const grid);

void mass_assignment_interlacing_cic(Catalogue const * const cat,
				     const Float x0[], const Float boxsize,
				     GridComplex* const grid);

#endif
