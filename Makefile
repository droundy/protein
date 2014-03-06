CXXFLAGS = -g -O2 -Wall -Werror

.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr

all: protein_microscopy

test-weights: test-weights.cpp weights.h


clean: rm -f protein_microscopy

protein_microscopy: protein_microscopy.cpp protein.h weights.h
