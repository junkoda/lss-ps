all: lib py

CC      := c++ -std=c++11
CXX     := c++ -std=c++11

OPT     := -DDOUBLEPRECISION
# OPT   += -fopenmp   # compiler option to enable OpenMP parallelisation

export CC CXX OPT

DIRS := lib py

.PHONY: lib py test clean check

lib:
	cd lib && $(MAKE) lib

py:
	cd py && $(MAKE) py

clean:
	for dir in $(DIRS); do (cd $$dir && $(MAKE) clean); done

