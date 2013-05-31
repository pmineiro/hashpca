CXXFLAGS=-Wall -O3 -funroll-loops -std=c++0x -march=native -mtune=native -ffast-math

all: pca

clean:
	rm -f pca

pca: pca.cc MurmurHash3.cpp $(wildcard *.hh)
	c++ $(CXXFLAGS) `pkg-config --cflags eigen3` $(word 1,$^) $(word 2,$^) -o $@
