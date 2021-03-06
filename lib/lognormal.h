#ifndef LOGRNORMAL_H
#define LOGRNORMAL_H 1

#include "grid.h"

void convert_P_to_delta_k(Grid* const grid,
			  const unsigned long seed,
			  const bool fix_amplitude=false);

#endif
