#!/bin/bash

# srun --job-name="time-arrow randst 18.60 28.60 94.00 full 120.00 1030.00" python pyplots/time-averaged-arrow.py randst 0.25 18.60 28.60 94.00 15.00 full_array 120.00 1030.00 &
# srun --job-name="time-arrow randst 18.60 28.60 94.00 full 220.00 1030.00" python pyplots/time-averaged-arrow.py randst 0.25 18.60 28.60 94.00 15.00 full_array 220.00 1030.00 &
# srun --job-name="time-arrow randst 18.60 28.60 94.00 full 320.00 1030.00" python pyplots/time-averaged-arrow.py randst 0.25 18.60 28.60 94.00 15.00 full_array 320.00 1030.00 &

# srun --job-name="time-arrow randst 18.60 28.60 94.00 exact 500.00 1030.00" python pyplots/time-averaged-arrow.py randst 0.25 18.60 28.60 94.00 15.00 exact 500.00 1030.00 &

# srun --job-name="time-aarrow randst 18.50 18.50 95.00 full_array 500.00 1100.00" python pyplots/time-averaged-arrow.py randst 0.25 18.50 18.50 95.00 15.00 full_array 500.00 1100.00 &
# srun --job-name="time-arrow randst 18.50 18.50 95.00 exact 500.00 1100.00" python pyplots/time-averaged-arrow.py randst 0.25 18.50 18.50 95.00 15.00 exact 500.00 1100.00 &

# srun --job-name="time-arrow randst 18.60 28.60 94.00 full_array 500.00 1100.00" python pyplots/time-averaged-arrow.py randst 0.25 18.60 28.60 94.00 15.00 full_array 500.00 1020.00 &
# srun --job-name="time-arrow randst 18.60 28.60 94.00 exact 500.00 1100.00" python pyplots/time-averaged-arrow.py randst 0.25 18.60 28.60 94.00 15.00 exact 500.00 1100.00 &

# srun --job-name="time-arrow stad 2.35 1.32 full 100.00 300.00" python pyplots/time-averaged-arrow.py stad 0.25 2.35 1.32 0.00 15.00 full_array 500.00 1200.00 &
# srun --job-name="time-arrow stad 2.35 1.32 exact 100.00 300.00" python pyplots/time-averaged-arrow.py stad 0.25 2.35 1.32 0.00 15.00 exact 500.00 1100.00 &
# srun --job-name="time-arrow stad 2.92 1.18 full 100.00 300.00" python pyplots/time-averaged-arrow.py stad 0.25 2.92 1.18 0.00 15.00 full_array 500.00 1200.00 &
# srun --job-name="time-arrow stad 2.92 1.18 exact 100.00 300.00" python pyplots/time-averaged-arrow.py stad 0.25 2.92 1.18 0.00 15.00 exact 500.00 1100.00 &

srun --job-name="time-arrow p 3.00 0.50 full 500.00 1100.00" python pyplots/time-averaged-arrow.py p 3.00 0.50 0.00 0.00 15.00 full_array 500.00 1100.00 &
srun --job-name="time-arrow p 3.00 0.50 exact 500.00 1100.00" python pyplots/time-averaged-arrow.py p 3.00 0.50 0.00 0.00 15.00 exact 500.00 1100.00 &


# Needed for main plot:
# 	data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/ave-time/contour-values-NflD-500-1900.dat \
# 	data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/ave-time/contour-values-NflD-500-980.dat \
# 	data/shape-stad/0_25-2_35-1_32-0_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-stad/0_25-2_35-1_32-0_00-15_00-full_array/ave-time/contour-values-NflD-500-1000.dat \
# 	data/shape-stad/0_25-2_92-1_18-0_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-stad/0_25-2_92-1_18-0_00-15_00-full_array/ave-time/contour-values-NflD-500-1000.dat \
# 	data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/ave-time/contour-values-NflD-500-1000.dat \
# 	data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-randst/0_25-18_50-18_50-95_00-15_00-exact/ave-time/contour-values-NflD-500-1800.dat \
# 	data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-randst/0_25-18_60-28_60-94_00-15_00-exact/ave-time/contour-values-NflD-500-1400.dat \
# 	data/shape-stad/0_25-2_35-1_32-0_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-stad/0_25-2_35-1_32-0_00-15_00-exact/ave-time/contour-values-NflD-500-780.dat \
# 	data/shape-stad/0_25-2_92-1_18-0_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-stad/0_25-2_92-1_18-0_00-15_00-exact/ave-time/contour-values-NflD-500-780.dat \
# 	data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/ave-time/ave-time-arrow-500-NflD.dat \
# 	data/shape-p/3_00-0_50-0_00-0_00-15_00-exact/ave-time/contour-values-NflD-500-1000.dat


# srun --job-name="time-arrow stad 3.00 1.18 full 500.00 1100.00" python pyplots/time-averaged-arrow.py stad 0.25 3.00 1.18 0.00 15.00 full_array 500.00 2000.00 &
# srun --job-name="time-arrow stad 3.00 1.18 exact 500.00 600.00" python pyplots/time-averaged-arrow.py stad 0.25 3.00 1.18 0.00 15.00 exact 500.00 2000.00 &


# srun --job-name="time-arrow p 3.50 0.50 exact 500.00 4000.00" python pyplots/time-averaged-arrow.py p 3.00 0.50 0.00 0.00 15.00 full_array 500.00 4000.00 &
# srun --job-name="time-arrow p 3.50 0.50 full_array 500.00 4000.00" python pyplots/time-averaged-arrow.py p 3.00 0.50 0.00 0.00 15.00 exact 500.00 4000.00 &
# srun --job-name="time-arrow p 3.00 0.50 exact 500.00 4000.00" python pyplots/time-averaged-arrow.py p 3.50 0.50 0.00 0.00 15.00 full_array 500.00 4000.00 &
# srun --job-name="time-arrow p 3.00 0.50 full_array 500.00 4000.00" python pyplots/time-averaged-arrow.py p 3.50 0.50 0.00 0.00 15.00 exact 500.00 4000.00

# srun --job-name="single-im-create p 3.00 0.50 full_array 300.00 500.00" python pyplots/single-image-creation.py p 3.00 0.50 0.00 0.00 15.00 full_array 300.00 500.00

