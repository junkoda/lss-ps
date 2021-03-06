#
# C++ library of large-scale structure power spectrum
#

default:
	cd .. && make

exec: lssps

LIB = liblssps.a liblssps.dylib
CPPFLAGS := -O3 $(OPENMP) $(OPT)


FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl
HDF5P_DIR ?= # parallel HDF5 library; e.g., brew install hdf5 --with-mpi

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR) $(HDF5P_DIR)

CPPFLAGS  += -I../lib $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += -L../lib $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

LIBS += -lgsl -lgslcblas -lboost_program_options -llssps 

OBJS := lssps.o
lssps: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@


#
# Dependeces
#
$(OBJS) : ../Makefile Makefile ../lib/liblssps.a

lssps.o: lssps.cpp ../lib/lssps.h ../lib/msg.h ../lib/config.h \
  ../lib/grid.h ../lib/catalogue.h ../lib/mass_assignment.h \
  ../lib/interlacing.h ../lib/power_spectrum.h

.PHONY: clean run dependence
clean:
	rm -f $(LIB) $(OBJS)
