#ifndef MEAN_DENSITY_H
#define MEAN_DENSITY_H 1

template<typename float_type>
void mean_density_from_grid(Grid const * const grid,
			    const Float fac,
			    const size_t np,
			    float_type const * xyz,
			    const size_t xyz_stride,
			    float_type* nbar,
			    const size_t nbar_stride)
{
  // Compute fac*grid at position xyz and set it to nbar array
  //
  // grid: Grid of scalar value
  // fac:  Multiplicative factor
  // np:   number of particles
  // xyz:  Array of positions
  // xyz_stride: bytes to next position
  // nbar: Output array
  // nbar_stide: bytes to next nbar data
  
  const int nc= grid->nc;
  const Float boxsize= grid->boxsize;
  const Float dx_inv= nc/boxsize;
  Float const * const x0 = grid->x0_box;
  //Float const * const d= grid->fx;

  //const size_t nckz= nc/2 + 1;
  
  for(size_t i=0; i<np; ++i) {
    int ix= (int) floor((xyz[0] - x0[0])*dx_inv);
    int iy= (int) floor((xyz[1] - x0[1])*dx_inv);
    int iz= (int) floor((xyz[2] - x0[2])*dx_inv);

    size_t iix= (ix + nc) % nc;
    size_t iiy= (iy + nc) % nc;
    size_t iiz= (iz + nc) % nc;

    //size_t index= (iix*nc + iiy)*nc + iiz;

    *nbar = fac*(*grid)(iix, iiy, iiz); //fac*d[index];
    
    xyz  = (float_type*) ((char*) xyz    + xyz_stride);
    nbar  = (float_type*) ((char*) nbar   + nbar_stride);
  }
}


#endif
