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
	paper/plot-ave.pdf \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/plots/box-plot_D.pdf \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/plots/box-plot_D.pdf \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/plots/box-plot_D.pdf \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/plots/box-plot_D.pdf \
	data/shape-stad/0_25-5_50-1_00-0_00-15_00-full_array/plots/box-plot_D.pdf \
	data/shape-stad/0_25-5_50-1_00-0_00-15_00-exact/plots/box-plot_D.pdf \
	data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/plots/single-image-plot.pdf \
	data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/plots/single-image-plot.pdf

paper/paper.pdf: paper/paper.tex ${ALL_FIGURES}
	echo ${ALL_FIGURES}
	cd paper && pdflatex paper.tex && bibtex paper && pdflatex paper.tex && pdflatex paper.tex

protein_test: protein_test.cpp protein_utils.cpp protein_membrane.cpp weights.cpp

# Below are rules to build each pdf file that needs building

paper/reactions.pdf: paper/reactions.svg
	inkscape --export-pdf $@ $<

# start time 29.501, period 45.002
data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/plots/single-image-plot.pdf: \
		pyplots/single-image-plot.py pyplots/mycolormap.py
	python $< p 3.00 0.50 0.00 0.00 15.00 full_array 312.00 350.00

data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/plots/single-image-plot.pdf: \
		pyplots/single-image-plot.py pyplots/mycolormap.py pyplots/single-image-creation.py \
		data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/nE/movie-frame-00600.dat
	python $< p 3.00 0.50 0.00 0.00 15.00 exact 266.00 304.00

data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/nE/movie-frame-00699.dat: pyplots/single-image-creation.py
	python pyplots/single-image-creation.py p 3.00 0.50 0.00 0.00 15.00 exact 312.00 350.00

paper/plot-ave.pdf: paper/plot-arrow-ave.py paper/mycolormap.py \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/ave-time/contour-values-NflD-500-850.dat \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/ave-time/contour-values-NflD-500-850.dat \
	data/shape-stad/0_25-2_35-1_32-0_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-stad/0_25-2_35-1_32-0_00-15_00-full_array/ave-time/contour-values-NflD-500-850.dat \
	data/shape-stad/0_25-2_92-1_18-0_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-stad/0_25-2_92-1_18-0_00-15_00-full_array/ave-time/contour-values-NflD-500-850.dat \
	data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/ave-time/contour-values-NflD-500-750.dat \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/ave-time/contour-values-NflD-500-850.dat \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/ave-time/contour-values-NflD-500-850.dat \
	data/shape-stad/0_25-2_35-1_32-0_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-stad/0_25-2_35-1_32-0_00-15_00-exact/ave-time/contour-values-NflD-500-850.dat \
	data/shape-stad/0_25-2_92-1_18-0_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-stad/0_25-2_92-1_18-0_00-15_00-exact/ave-time/contour-values-NflD-500-850.dat \
	data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
	data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/ave-time/contour-values-NflD-500-750.dat
	python paper/plot-arrow-ave.py

data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/plots/box-plot_D.pdf: pyplots/box_plot.py \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/box-plot.dat \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/sections.dat
	python pyplots/box_plot.py randst 0.25 18.50 18.50 95.00 15.00 exact 500 1000

data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/plots/box-plot_D.pdf: pyplots/box_plot.py \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/box-plot.dat \
	data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/sections.dat
	python pyplots/box_plot.py randst 0.25 18.50 18.50 95.00 15.00 full_array 500 1000

data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/plots/box-plot_D.pdf: pyplots/box_plot.py \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/box-plot.dat \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/sections.dat
	python pyplots/box_plot.py randst 0.25 18.60 28.60 94.00 15.00 exact 500 1000

data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/plots/box-plot_D.pdf: pyplots/box_plot.py \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/box-plot.dat \
	data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/sections.dat
	python pyplots/box_plot.py randst 0.25 18.60 28.60 94.00 15.00 full_array 500 1000

data/shape-stad/0_25-5_50-1_00-0_00-15_00-exact/plots/box-plot_D.pdf: pyplots/box_plot.py \
	data/shape-stad/0_25-5_50-1_00-0_00-15_00-exact/box-plot.dat \
	data/shape-stad/0_25-5_50-1_00-0_00-15_00-exact/sections.dat
	python pyplots/box_plot.py stad 0.25 5.50 1.00 0.00 15.00 exact 1500 2100

data/shape-stad/0_25-5_50-1_00-0_00-15_00-full_array/plots/box-plot_D.pdf: pyplots/box_plot.py \
	data/shape-stad/0_25-5_50-1_00-0_00-15_00-full_array/box-plot.dat \
	data/shape-stad/0_25-5_50-1_00-0_00-15_00-full_array/sections.dat
	python pyplots/box_plot.py stad 0.25 5.50 1.00 0.00 15.00 full_array 1500 2100
