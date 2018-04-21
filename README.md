# lss-ps

A flexible toolkit to compute power spectrum for large-scale structure
cosmology with Python


## Setup

### Requirements

* C++11 compiler
* Python3 and Numpy
* FFTW and GSL libraries

Edit Makefile for C++ compiler, non-standard library path.

```bash
$ make
```

This compiles the code and installs the library.

If you do not have root priviledge to install Python3 or the python
library you can try [miniconda](https://conda.io/miniconda.html).

### Check

To check if `lssps` is installed, import it from Python3

```bash
$ python3
>>> import lssps
```


## Example

```python

import numpy as np
import lssps

# load 4 columns x, y, z, nbar
# where nbar is the mean density (density without clustering) at each position
data = np.loadtxt('data.txt') # load data xyz coordinate
rand = np.loadtxt('rand.txt') # load random catalogue xyz coordinate

# Assign density to grids
boxsize = 3400.0       # length of a cube on a side
x0 = (0.0, 0.0, 0.0)   # xyz coordinate of the corner of the box
                       # [x0, x0 + boxsize] is expected to contain all the
		       # objects

grid = lssps.grid.zeros(nc, boxsize, x0, offset=0.5, interlacing=True)
grid_rand = lssps.grid.zeros(nc, boxsize, x0, offset=0.5, interlacing=True)

grid.assign_density(xyz=data[:, 0:3], nbar=data[:, 3], mas='CIC')
grid_rand.assign_density(xyz=rand[:, 0:3], nbar=rand[:, 3], mas='CIC')
# you can also set weight=...

# Compute fluctuation field F(x) = n_data(x) - alpha*n_rand(x)
# where alpha = data/random ratio
grid.compute_fluctuation(grid_rand)

# compute_fluctuation can be called without randoms for homogeneous data in
# a periodic box (e.g. simulation snapthot)
# grid.compute_fluctuation()


# Compute power spectrum
ps = lssps.compute_plane_parallel(grid, k_min=0.0, k_max=1.0, dk=0.01)

#
# You can also use FFT based Yamamoto estimators for wide-angle survey
#
# ps = lssps.power_spectrum.compute_yamamoto(grid, 'scoccimarro',
#                                            k_min=0.0, k_max=1.0, dk=0.01)
# or,
# ps = lssps.power_spectrum.compute_yamamoto(grid, 'bianchi',
#                                            k_min=0.0, k_max=1.0, dk=0.01)
#

# Print power spectrum
# quadrupole and hexadecapole, ps.P2, ps.P4, are also available
for i in range(len(ps)):
    if ps.nmodes[i] > 0:
         print('%e %e %d' % (ps.k[i], ps.P0[i], ps.nmodes[i]))
```
