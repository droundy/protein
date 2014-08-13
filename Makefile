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
	data/shape-randst/plots/box-plot_D--randst-25-800-800-9500-1500-exact.pdf \
	data/shape-randst/plots/box-plot_D--randst-25-1000-1700-9400-1500-exact.pdf \
	data/shape-stad/plots/plot-time-averaged-arrow-500-520-NflD-stad-25-400-150-0-1500-full_array.pdf \
	data/shape-p/plots/image-plot--p-300-50-0-0-1500-exact.pdf \
	data/shape-p/plots/single-image-plot--p-300-50-0-0-1500-full_array.pdf \
	data/shape-stad/plots/time-map-NflD-stad-25-600-200-0-1500-full_array.pdf \
	data/shape-stad/plots/time-map-NflD-stad-25-400-150-0-1500-full_array.pdf \
	data/shape-randst/plots/time-map-NflD-randst-25-1850-1850-9500-1500-full_array.pdf \
	data/shape-randst/plots/time-map-NflD-randst-25-2500-2500-9500-1500-full_array.pdf \
	data/shape-randst/plots/time-map-NflD-randst-25-1860-2860-9400-1500-full_array.pdf \
	data/shape-randst/plots/time-map-NflD-randst-25-2300-3200-9400-1500-full_array.pdf \
	paper/plot-ave.pdf

paper/paper.pdf: paper/paper.tex ${ALL_FIGURES}
	echo ${ALL_FIGURES}
	cd paper && pdflatex paper.tex && bibtex paper && pdflatex paper.tex && pdflatex paper.tex

protein_test: protein_test.cpp protein_utils.cpp protein_membrane.cpp weights.cpp

# Below are rules to build each pdf file that needs building

paper/reactions.pdf: paper/reactions.svg
	inkscape --export-pdf $@ $<

# start time 29.501, period 45.002
data/shape-p/plots/image-plot--p-300-50-0-0-1500-exact.pdf: pyplots/image_plot.py
	mkdir -p data/shape-p/plots
	python $< p 3.00 0.50 0.00 0.00 15.00 exact 266.00 304.00

data/shape-p/plots/single-image-plot--p-300-50-0-0-1500-full_array.pdf: pyplots/single-image-plot.py
	python $< p 3.00 0.50 0.00 0.00 15.00 full_array 330.00 368.00

paper/plot-ave.pdf: paper/plot-arrow-ave.py \
	./data/shape-randst/plots/ave-time/maxima-arrowNflD-randst-0.25-18.50-18.50-95.00-15.00-full_array.dat \
	./data/shape-randst/plots/ave-time/contour-values-300-1500-NflD-randst-0.25-18.50-18.50-95.00-15.00-full_array.dat
	python paper/plot-arrow-ave.py

data/shape-stad/plots/plot-time-averaged-arrow-500-520-NflD-stad-25-400-150-0-1500-full_array.pdf: pyplots/plot-ave-arrow.py \
	./data/shape-stad/plots/ave-time/maxima-arrow-500-NflD-stad-0.25-4.00-1.50-0.00-15.00-full_array.dat \
	./data/shape-stad/plots/ave-time/contour-values-500-520-NflD-stad-0.25-4.00-1.50-0.00-15.00-full_array.dat
	python pyplots/plot-ave-arrow.py stad 0.25 4.00 1.50 0.00 15.00 full_array 500 520

data/shape-randst/plots/box-plot_D--randst-25-800-800-9500-1500-exact.pdf: pyplots/box_plot.py data/shape-randst/box-plot--randst-0.25-8.00-8.00-95.00-15.00.dat
	python pyplots/box_plot.py randst 0.25 8.00 8.00 95.00 15.00 exact 0 5

data/shape-randst/plots/box-plot_D--randst-25-1000-1700-9400-1500-exact.pdf: pyplots/box_plot.py data/shape-randst/box-plot--randst-0.25-10.00-17.00-94.00-15.00.dat
	python pyplots/box_plot.py randst 0.25 10.00 17.00 94.00 15.00 exact 0 5

data/shape-stad/plots/time-map-NflD-stad-25-600-200-0-1500-full_array.pdf: pyplots/time_map.py \
	data/shape-stad/time-map-NflD-stad-0.25-6.00-2.00-0.00-15.00-full_array.dat
	python pyplots/time_map.py stad 0.25 6.00 2.00 0.00 15.00 full_array 0.00 10.00
data/shape-stad/plots/time-map-NflD-stad-25-400-150-0-1500-full_array.pdf: pyplots/time_map.py \
	data/shape-stad/time-map-NflD-stad-0.25-4.00-1.50-0.00-15.00-full_array.dat
	python pyplots/time_map.py stad 0.25 4.00 1.50 0.00 15.00 full_array 0.00 10.00
data/shape-randst/plots/time-map-NflD-randst-25-1850-1850-9500-1500-full_array.pdf: pyplots/time_map.py \
	data/shape-randst/time-map-NflD-randst-0.25-18.50-18.50-95.00-15.00-full_array.dat
	python pyplots/time_map.py randst 0.25 18.50 18.50 95.00 15.00 full_array 0.00 10.00
data/shape-randst/plots/time-map-NflD-randst-25-2500-2500-9500-1500-full_array.pdf: pyplots/time_map.py \
	data/shape-randst/time-map-NflD-randst-0.25-25.00-25.00-95.00-15.00-full_array.dat
	python pyplots/time_map.py randst 0.25 25.00 25.00 95.00 15.00 full_array 0.00 10.00
data/shape-randst/plots/time-map-NflD-randst-25-1860-2860-9400-1500-full_array.pdf: pyplots/time_map.py \
	data/shape-randst/time-map-NflD-randst-0.25-18.60-28.60-94.00-15.00-full_array.dat
	python pyplots/time_map.py randst 0.25 18.60 28.60 94.00 15.00 full_array 0.00 10.00
data/shape-randst/plots/time-map-NflD-randst-25-2300-3200-9400-1500-full_array.pdf: pyplots/time_map.py \
	data/shape-randst/time-map-NflD-randst-0.25-23.00-32.00-94.00-15.00-full_array.dat
	python pyplots/time_map.py randst 0.25 23.00 32.00 94.00 15.00 full_array 0.00 10.00


#data/shape-p/plots/plot-time-averaged-arrow-300-980-NflD-p-400-50-0-0-1500-full_array.pdf
#data/shape-randst/plots/plot-time-averaged-arrow-300-1500-NflD-randst-25-1850-1850-9500-1500-full_array.pdf
#data/shape-stad/plots/plot-time-averaged-arrow-500-1500-NflD-stad-25-400-150-0-1500-full_array.pdf
