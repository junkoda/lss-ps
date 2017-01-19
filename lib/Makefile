#
# C++ library of large-scale structure power spectrum
#

default:
	cd .. && make lib

LIB = liblssps.a liblssps.dylib
CPPFLAGS := -O3 $(OPENMP) $(OPT)
#OPT     := -DDOUBLEPRECISION # debug!!!

lib: $(LIB)

FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl
HDF5P_DIR ?= # parallel HDF5 library; e.g., brew install hdf5 --with-mpi

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR) $(HDF5P_DIR)

CPPFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

LIBS += -lgsl -lgslcblas

OBJS := msg.o config.o catalogue.o grid.o mass_assignment.o power_spectrum.o
OBJS += interlacing.o


ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  # Single precision FFTW
  FFTWSUF=f
endif
LIBS += -lfftw3$(FFTWSUF)

ifdef OPENMP
  LIBS += -lfftw3$(FFTWSUF)_omp
  #LIBS += -lfftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif

#
# Dependeces
#
$(OBJS) : ../Makefile Makefile

#
# Library libfs.a
#
liblssps.a: $(OBJS)
	ar r $@ $(OBJS)

liblssps.dylib: $(OBJS)
	echo "OPT= $(OPT)"
	$(CXX) -dynamiclib $(LIBS) $(OBJS) -o $@

.PHONY: clean run dependence
clean:
	rm -f $(LIB) $(OBJS)
