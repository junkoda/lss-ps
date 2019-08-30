default: py

# Set C++11 compiler
#
# e.g., CXX = g++ -fopenmp -std=c++11
#
# 1. Uses $CXXOPENMP variable if it is defined;
# 2. CXX = $(CXXOPENMP) $(OPENMP)-std=c++11 
# 3. if CXXOPENMP is not defined, uses g++
# Set OpenMP flag to enable parallelisation

OPENMP  := #-fopenmp

CXXOPENMP = g++
CXX     := $(CXXOPENMP) $(OPENMP) -std=c++11
CC      := $(CXX)

# Extra compile options
# Use double pricision; comment out to use single precision (float in C)
OPT     := -DDOUBLEPRECISION 


# Set library path if they are in non-standard directories
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR)

IDIRS    += $(foreach dir, $(DIR_PATH), $(dir)/include)
LDIRS    += $(foreach dir, $(DIR_PATH), $(dir)/lib)

LIBS += gsl gslcblas

# A single/double precision variable for FFTW library
ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  # Single precision FFTW
  FFTWSUF=f
endif
LIBS += fftw3$(FFTWSUF)

ifdef OPENMP
  # FFTW can be parallelised with OpenMP or threads
  LIBS += fftw3$(FFTWSUF)_omp
  #LIBS += fftw3$(FFTWSUF)_threads # replace the line above for thread parallelization instead of omp
endif

# These evironment variables are used in py/setup.py
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
