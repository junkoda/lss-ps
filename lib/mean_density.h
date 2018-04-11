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

  for(size_t i=0; i<np; ++i) {
    /*
    int ix= (int) floor((xyz[0] - x0[0])*dx_inv);
    int iy= (int) floor((xyz[1] - x0[1])*dx_inv);
    int iz= (int) floor((xyz[2] - x0[2])*dx_inv);

    size_t iix= (ix + nc) % nc;
    size_t iiy= (iy + nc) % nc;
    size_t iiz= (iz + nc) % nc;

    //size_t index= (iix*nc + iiy)*nc + iiz;

    
    

    *nbar = fac*(*grid)(iix, iiy, iiz); //fac*d[index];
    */

    // TSC interpolation
    int ix[3], ix0[3], ix1[3], ix2[3];
    double w0[3], w1[3], w2[3];

    for(int k=0; k<3; ++k) {
      Float rx = (xyz[k] - x0[k])*dx_inv;
      ix[k] = (int) floor(rx + 0.5);
      ix0[k]= (ix[k] - 1 + grid->nc) % grid->nc; // left grid point (periodic)
      ix1[k]= (ix[k]     + grid->nc) % grid->nc; // nearest grid point (periodic)
      ix2[k]= (ix[k] + 1 + grid->nc) % grid->nc; // right grid point (periodic)

      double dx1 = rx - ix[k];
      double dx2 = 0.5 - dx1;
      
      w0[k] = 0.5*dx2*dx2;
      w1[k] = 0.75 - dx1*dx1;
      w2[k] = 0.25 - 0.5*dx2*dx2 + dx1*dx1;
    }

    double n= 0.0;

    n += (*grid)(ix0[0], ix0[1], ix0[2])*w0[0]*w0[1]*w0[2];
    n += (*grid)(ix0[0], ix0[1], ix1[2])*w0[0]*w0[1]*w1[2];
    n += (*grid)(ix0[0], ix0[1], ix2[2])*w0[0]*w0[1]*w2[2];
    n += (*grid)(ix0[0], ix1[1], ix0[2])*w0[0]*w1[1]*w0[2];
    n += (*grid)(ix0[0], ix1[1], ix1[2])*w0[0]*w1[1]*w1[2];
    n += (*grid)(ix0[0], ix1[1], ix2[2])*w0[0]*w1[1]*w2[2];
    n += (*grid)(ix0[0], ix2[1], ix0[2])*w0[0]*w2[1]*w0[2];
    n += (*grid)(ix0[0], ix2[1], ix1[2])*w0[0]*w2[1]*w1[2];
    n += (*grid)(ix0[0], ix2[1], ix2[2])*w0[0]*w2[1]*w2[2];

    n += (*grid)(ix1[0], ix0[1], ix0[2])*w1[0]*w0[1]*w0[2];
    n += (*grid)(ix1[0], ix0[1], ix1[2])*w1[0]*w0[1]*w1[2];
    n += (*grid)(ix1[0], ix0[1], ix2[2])*w1[0]*w0[1]*w2[2];
    n += (*grid)(ix1[0], ix1[1], ix0[2])*w1[0]*w1[1]*w0[2];
    n += (*grid)(ix1[0], ix1[1], ix1[2])*w1[0]*w1[1]*w1[2];
    n += (*grid)(ix1[0], ix1[1], ix2[2])*w1[0]*w1[1]*w2[2];
    n += (*grid)(ix1[0], ix2[1], ix0[2])*w1[0]*w2[1]*w0[2];
    n += (*grid)(ix1[0], ix2[1], ix1[2])*w1[0]*w2[1]*w1[2];
    n += (*grid)(ix1[0], ix2[1], ix2[2])*w1[0]*w2[1]*w2[2];

    n += (*grid)(ix2[0], ix0[1], ix0[2])*w2[0]*w0[1]*w0[2];
    n += (*grid)(ix2[0], ix0[1], ix1[2])*w2[0]*w0[1]*w1[2];
    n += (*grid)(ix2[0], ix0[1], ix2[2])*w2[0]*w0[1]*w2[2];
    n += (*grid)(ix2[0], ix1[1], ix0[2])*w2[0]*w1[1]*w0[2];
    n += (*grid)(ix2[0], ix1[1], ix1[2])*w2[0]*w1[1]*w1[2];
    n += (*grid)(ix2[0], ix1[1], ix2[2])*w2[0]*w1[1]*w2[2];
    n += (*grid)(ix2[0], ix2[1], ix0[2])*w2[0]*w2[1]*w0[2];
    n += (*grid)(ix2[0], ix2[1], ix1[2])*w2[0]*w2[1]*w1[2];
    n += (*grid)(ix2[0], ix2[1], ix2[2])*w2[0]*w2[1]*w2[2];

    *nbar = fac*n;

    xyz  = (float_type*) ((char*) xyz    + xyz_stride);
    nbar  = (float_type*) ((char*) nbar   + nbar_stride);
  }
}


#endif
