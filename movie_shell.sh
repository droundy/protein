#!/bin/bash

################################################
#From that paper plot, to get the times right:
# arg_set = ["randst/0.25-18.50-18.50-95.00-15.00-exact",
#            "randst/0.25-18.50-18.50-95.00-15.00-full_array",
#            "randst/0.25-18.60-28.60-94.00-15.00-exact",
#            "randst/0.25-18.60-28.60-94.00-15.00-full_array",
#            "stad/0.25-2.35-1.32-0.00-15.00-exact",
#            "stad/0.25-2.35-1.32-0.00-15.00-full_array",
#            "stad/0.25-2.92-1.18-0.00-15.00-exact",
#            "stad/0.25-2.92-1.18-0.00-15.00-full_array",
#            "p/3.00-0.50-0.00-0.00-15.00-exact",
#            "p/3.00-0.50-0.00-0.00-15.00-full_array"]
#
# bound_times = [500,850,500,850,500,850,500,850,500,850,500,850,500,850,500,850,500,750,500,750]
#
###################################################


echo starting!
# srun python pyplots/movie-image-creation.py p 1.00 0.50 0.00 0.00 15.00 full_array 300.00 1000.00 &
# srun python pyplots/movie-image-creation.py p 2.00 0.50 0.00 0.00 15.00 full_array 300.00 1000.00 &

# data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/plots/correlation.pdf: pyplots/correlation-plot.py \
# data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/plots/correlation.pdf: pyplots/correlation-plot.py \
# data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/plots/correlation.pdf: pyplots/correlation-plot.py \
# data/shape-stad/0_25-2_92-1_18-0_00-15_00-full_array/plots/correlation.pdf: pyplots/correlation-plot.py \
# data/shape-stad/0_25-2_35-1_32-0_00-15_00-full_array/plots/correlation.pdf: pyplots/correlation-plot.py \


# srun --job-name="mov-im-create p 3.00 0.50 full 500 1000" python pyplots/movie-image-creation.py p 3.00 0.50 0.00 0.00 15.00 full_array 500.00 750.00 &
# srun --job-name="mov-im-create randst 18.60 28.60 full 500 1000" python pyplots/movie-image-creation.py randst 0.40 18.60 28.60 94.00 15.00 full_array 500.00 850.00 &
# srun --job-name="mov-im-create randst 18.50 18.50 full 500 1000" python pyplots/movie-image-creation.py randst 0.40 18.50 18.50 95.00 15.00 full_array 500.00 850.00 &
# srun --job-name="mov-im-create stad 2.92 1.18 full 500 1000" python pyplots/movie-image-creation.py stad 0.40 2.92 1.18 0.00 15.00 full_array 500.00 850.00 &
# srun --job-name="mov-im-create stad 2.35 1.32 full 500 1000" python pyplots/movie-image-creation.py stad 0.40 2.35 1.32 0.00 15.00 full_array 500.00 850.00 &

# srun --job-name="mov-maker p 3.00 0.50 full 500 750" python pyplots/movie-maker.py p 3.00 0.50 0.00 0.00 15.00 full_array 500.00 750.00 &
# srun --job-name="mov-maker randst 18.60 28.60 full 500 850" python pyplots/movie-maker.py randst 0.40 18.60 28.60 94.00 15.00 full_array 500.00 850.00 &
# srun --job-name="mov-maker randst 18.50 18.50 full 500 850" python pyplots/movie-maker.py randst 0.40 18.50 18.50 95.00 15.00 full_array 500.00 850.00 &
# srun --job-name="mov-maker stad 2.92 1.18 full 500 850" python pyplots/movie-maker.py stad 0.40 2.92 1.18 0.00 15.00 full_array 500.00 850.00 &
# srun --job-name="mov-maker stad 2.35 1.32 full 500 850" python pyplots/movie-maker.py stad 0.40 2.35 1.32 0.00 15.00 full_array 500.00 850.00 &
