SHELL=/bin/zsh

CXXFLAGS=-O3 -funroll-loops -std=c++0x -march=native -mtune=native -ffast-math
EIGENCFLAGS := $(shell pkg-config --cflags eigen3)

all: pca

clean:
	rm -f pca

pca: pca.cc MurmurHash3.cpp $(wildcard *.hh)
	c++ $(CXXFLAGS) $(EIGENCFLAGS) $(word 1,$^) $(word 2,$^) -o $@

check: testsvd.ok testpca.ok

%.ok: %
	@printf "%s" "=== running test $(word 1,$^) ===" 1>&2
	@./$(word 1,$^)
	@echo "=== $(word 1,$^) passed ===" 1>&2
