CC      := c++ -std=c++11
CXX     := c++ -std=c++11

OPT     := -DDOUBLEPRECISION

export CC CXX OPT

.PHONY: lib py test clean check

lib:
	cd lib && $(MAKE) lib

py:
	cd py && $(MAKE) py

clean:
	for dir in $(DIRS); do (cd $$dir && $(MAKE) clean); done

