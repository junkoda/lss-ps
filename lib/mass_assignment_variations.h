#ifndef MASS_ASSIGNMENT_VARIATIONS_H
#define MASS_ASSIGNMENT_VARIATIONS_H 1

void mass_assignment_from_array_serial(double const * const xyz,
				       const size_t xyz_stride,
				       double const * const weight,
				       const size_t weight_stride,
				       double const * const nbar,
				       const size_t nbar_stride,
				       const size_t np,
				       const int mas,
				       bool parallelise,
				       Grid* const grid);

void mass_assignment_from_array_atomic(double const * const xyz,
				       const size_t xyz_stride,
				       double const * const weight,
				       const size_t weight_stride,
				       double const * const nbar,
				       const size_t nbar_stride,
				       const size_t np,
				       const int mas,
				       bool parallelise,
				       Grid* const grid);

#endif
