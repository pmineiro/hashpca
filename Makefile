VEEDUBPARSE=../veedubparse
CXXFLAGS=-Wall -O3 -funroll-loops -std=c++0x -march=native -mtune=native -ffast-math

all: veedubparsecheck pca

clean:
	rm -f pca

check: veedubparsecheck pca
	cd tests && $(MAKE) check

veedubparsecheck:
	@test -f $(VEEDUBPARSE)/MurmurHash3.cpp || {				\
	  echo "ERROR: you need to tell me where to find veedubparse" 1>&2;	\
	  echo "       ( available from https://github.com/pmineiro/veedubparse )" 1>&2; \
	  echo "       I'm currently looking for it at: $(VEEDUBPARSE)" 1>&2; \
	  echo "       edit the Makefile and try again" 1>&2; \
	}

pca: pca.cc $(VEEDUBPARSE)/MurmurHash3.cpp $(wildcard *.hh) 
	c++ $(CXXFLAGS) `pkg-config --cflags eigen3` -I $(VEEDUBPARSE) $(word 1,$^) $(word 2,$^) -o $@
