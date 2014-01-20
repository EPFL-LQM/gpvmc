CXX=g++
CFLAGS=-g3 -O2 -Wall -std=c++0x

#all: concatenate
all: utils concatenate

utils: vmc.py stagflux.py setup.py
	python setup.py sdist

concatenate: concatenate.o
	$(CXX) -o concatenate concatenate.o -lhdf5 -lhdf5_hl

concatenate.o: concatenate.cpp
	$(CXX) -c -o concatenate.o concatenate.cpp $(CFLAGS)

install: utils
	python setup.py install --user
