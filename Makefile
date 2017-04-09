all: lib py #exec

CC      := c++ -std=c++11
CXX     := c++ -std=c++11

OPT     := -DDOUBLEPRECISION
# OPT   += -fopenmp   # compiler option to enable OpenMP parallelisation

export CC CXX OPT

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
