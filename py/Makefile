# make
#   compile this python module
# make install
#   install the module to the python module directory

default:
	cd .. && $(MAKE) py

py: build install

#CC   ?= mpic++ -std=c++11
#SRC  := py_package.cpp py_msg.cpp
#CFLAGS = -DABC

build: $(SRC)
	python3 setup.py config build_ext --inplace

install:
	python3 setup.py install


.PHONY: clean check

clean:
	rm -rf build fs/*.so

check:
	cd fs && $(MAKE) check

