CXXFLAGS = -g -O2 -Wall -Werror

.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr

all: protein_microscopy

test-weights: test-weights.cpp weights.h weights.cpp

clean: rm -f protein_microscopy protein_test

protein_microscopy: protein_microscopy.cpp protein_weights.cpp protein_utils.cpp \
	 protein_membrane.cpp weights.cpp

ALL_FIGURES = \
	data/shape-p/plots/box-plot_D--p-300-50-0-0-1500.pdf

paper/paper.pdf: paper/paper.tex \
		paper/reactions.pdf data/shape-p/plots/image-plot--p-300-50-0-0-1500.pdf \
		data/shape-randst/plots/paper-arrow-plot.pdf \
		${ALL_FIGURES}
	echo ${ALL_FIGURES}
	cd paper && pdflatex paper.tex && bibtex paper && pdflatex paper.tex && pdflatex paper.tex

paper/reactions.pdf: paper/reactions.svg
	inkscape --export-pdf $@ $<

${ALL_FIGURES}: batch jobs $(wildcard data/shape-*/*.dat) $(wildcard pyplots/*.py)
	./batch plot include="-paper"


# start time 29.501, period 45.002
data/shape-p/plots/image-plot--p-300-50-0-0-1500.pdf: pyplots/image_plot.py
	mkdir -p data/shape-p/plots
	python $< p 3.00 0.50 0.00 0.00 15.00 266.00 304.00

data/shape-randst/plots/paper-arrow-plot.pdf: paper/arrow-plot.py $(wildcard data/shape-randst/*.dat)
	mkdir -p data/shape-randst/plots
	python $<

protein_test: protein_test.cpp protein_utils.cpp protein_membrane.cpp weights.cpp

