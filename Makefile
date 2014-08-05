CXXFLAGS = -g -O2 -Wall -Werror

.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr

all: protein_microscopy paper/paper.pdf test-weights

test-weights: test-weights.cpp weights.cpp weights.h
	g++ ${CXXFLAGS} -o test-weights test-weights.cpp weights.cpp

clean: rm -f protein_microscopy protein_test

protein_microscopy: protein_microscopy.cpp protein_weights.cpp protein_utils.cpp \
	 protein_membrane.cpp weights.cpp

# ALL_FIGURES should be a list of all the pdf files that are required
# to build our paper.  Or jpg files or whatever.
ALL_FIGURES = \
	paper/reactions.pdf \
	data/shape-randst/plots/box-plot_D--randst-25-800-800-9500-1500.pdf \
	data/shape-p/plots/image-plot--p-300-50-0-0-1500.pdf

paper/paper.pdf: paper/paper.tex ${ALL_FIGURES}
	echo ${ALL_FIGURES}
	cd paper && pdflatex paper.tex && bibtex paper && pdflatex paper.tex && pdflatex paper.tex

protein_test: protein_test.cpp protein_utils.cpp protein_membrane.cpp weights.cpp

# Below are rules to build each pdf file that needs building

paper/reactions.pdf: paper/reactions.svg
	inkscape --export-pdf $@ $<

# start time 29.501, period 45.002
data/shape-p/plots/image-plot--p-300-50-0-0-1500.pdf: pyplots/image_plot.py fixme-add-data-here.dat
	mkdir -p data/shape-p/plots
	python $< p 3.00 0.50 0.00 0.00 15.00 266.00 304.00

data/shape-randst/plots/paper-arrow-plot.pdf: paper/arrow-plot.py fixme-add-data-here.dat
	mkdir -p data/shape-randst/plots
	python $<

data/shape-randst/plots/box-plot_D--randst-25-800-800-9500-1500.pdf: pyplots/box_plot.py fix-me-add-files-here.dat
	python pyplots/box_plot.py randst .25 8.00 8.00 95.00 15.00 exact 0 10
