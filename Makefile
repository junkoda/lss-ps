default: py

# Set C++11 compiler
CXX     := c++ -std=c++11
CC      := $(CXX)

OPT     := -DDOUBLEPRECISION
# OPT   += -fopenmp   # compiler option to enable OpenMP parallelisation

# Set library path if they are in non-standard directories
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR)

IDIRS    += $(foreach dir, $(DIR_PATH), $(dir)/include)
LDIRS    += $(foreach dir, $(DIR_PATH), $(dir)/lib)

LIBS += gsl gslcblas

ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  # Single precision FFTW
  FFTWSUF=f
endif
LIBS += fftw3$(FFTWSUF)

ifdef OPENMP
  LIBS += fftw3$(FFTWSUF)_omp
  #LIBS += fftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif


export CC CXX OPT IDIRS LDIRS LIBS 


DIRS := lib py

.PHONY: lib py exec test clean check doc

lib:
	cd lib && $(MAKE) lib

py:
	cd py && $(MAKE) py

exec:
	cd exec && $(MAKE) exec

clean:
	for dir in $(DIRS); do (cd $$dir && $(MAKE) clean); done

doc:
	cd doc && $(MAKE) html
	@echo "Documents generated in doc/_build/html"
