# lss-ps
Power spectrum estimator for large-scale structure cosmology



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

```bash
$ python3
>>> import lssps
```


