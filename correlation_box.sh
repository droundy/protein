#!/bin/bash

echo starting !

srun --job-name="corr randst 18.50 18.50 95.00 full_array 0.00 4500.00" python pyplots/density-correlation.py randst 0.25 18.50 18.50 95.00 15.00 full_array 0.00 &
srun --job-name="corr randst 18.50 18.50 95.00 exact 0.00 5000.00" python pyplots/density-correlation.py randst 0.25 18.50 18.50 95.00 15.00 exact 0.00 &

srun --job-name="corr randst 18.60 28.60 94.00 full_array 0.00 1020.00" python pyplots/density-correlation.py randst 0.25 18.60 28.60 94.00 15.00 full_array 0.00 &
srun --job-name="corr randst 18.60 28.60 94.00 exact 0.00 1800.00" python pyplots/density-correlation.py randst 0.25 18.60 28.60 94.00 15.00 exact 0.00 &

srun --job-name="corr stad 2.35 1.32 full 0.00 10000.00" python pyplots/density-correlation.py stad 0.25 2.35 1.32 0.00 15.00 full_array 0.00 &
srun --job-name="corr stad 2.35 1.32 exact 0.00 6000.00" python pyplots/density-correlation.py stad 0.25 2.35 1.32 0.00 15.00 exact 0.00 &

srun --job-name="corr stad 2.92 1.18 full 0.00 1410.00" python pyplots/density-correlation.py stad 0.25 2.92 1.18 0.00 15.00 full_array 0.00 &
srun --job-name="corr stad 2.92 1.18 exact 0.00 910.00" python pyplots/density-correlation.py stad 0.25 2.92 1.18 0.00 15.00 exact 0.00 &

srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 &
srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00

####################################################
# srun --job-name="corr randst 18.50 18.50 95.00 full_array 0.00 10.00" python pyplots/correlation-plot.py randst 0.25 18.50 18.50 95.00 15.00 full_array 0.00 1125.00 rl &
# srun --job-name="corr randst 18.50 18.50 95.00 exact 0.00 10.00" python pyplots/correlation-plot.py randst 0.25 18.50 18.50 95.00 15.00 exact 0.00 1125.00 rl &

# srun --job-name="corr randst 18.60 28.60 94.00 full_array 0.00 10.00" python pyplots/correlation-plot.py randst 0.25 18.60 28.60 94.00 15.00 full_array 0.00 255.00 rl &
# srun --job-name="corr randst 18.60 28.60 94.00 exact 0.00 10.00" python pyplots/correlation-plot.py randst 0.25 18.60 28.60 94.00 15.00 exact 0.00 450.00 rl &

# srun --job-name="corr stad 2.35 1.32 full 0.00 10.00" python pyplots/correlation-plot.py stad 0.25 2.35 1.32 0.00 15.00 full_array 0.00 1500.00 rl &
# srun --job-name="corr stad 2.35 1.32 exact 0.00 10.00" python pyplots/correlation-plot.py stad 0.25 2.35 1.32 0.00 15.00 exact 0.00 1500.00 rl &

# srun --job-name="corr stad 2.92 1.18 full 0.00 10.00" python pyplots/correlation-plot.py stad 0.25 2.92 1.18 0.00 15.00 full_array 0.00 225.00 rl &
# srun --job-name="corr stad 2.92 1.18 exact 0.00 10.00" python pyplots/correlation-plot.py stad 0.25 2.92 1.18 0.00 15.00 exact 0.00 225.00 rl &

# srun --job-name="corr p 3.00 0.50 full 0.00 10.00" python pyplots/correlation-plot.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 1500.00 rl &
# srun --job-name="corr p 3.00 0.50 exact 0.00 10.00" python pyplots/correlation-plot.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 1500.00 rl


########################

# list=(yes no)
# for i in $list;
# do
#     echo $i
# done

# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 600.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 1000.00 &
# cd ../new-protein-0
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 300.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 1000.00 &
# cd ../new-protein-1
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 1000.00 &
# cd ../new-protein-2
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &
# cd ../new-protein-3
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &
# cd ../new-protein-4
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &
# cd ../new-protein-5
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &
# cd ../new-protein-6
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &
# cd new-protein-7
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &
# cd new-protein-8
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &
# cd new-protein-9
# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 200.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 200.00 &

echo finished!!!


# python pyplots/correlation-plot.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 400.00 rl 5
# python pyplots/correlation-plot.py randst 0.25 18.60 28.60 94.00 15.00 full_array 0.00 220.00 rl 2
# python pyplots/correlation-plot.py randst 0.25 18.50 18.50 95.00 15.00 full_array 0.00 450.00 rl 3
# python pyplots/correlation-plot.py stad 0.25 2.92 1.18 0.00 15.00 full_array 0.00 200.00 rl 2
# python pyplots/correlation-plot.py stad 0.25 2.35 1.32 0.00 15.00 full_array 0.00 380.00 rl 5

# data/shape-p/3_00-0_50-0_00-0_00-15_00-full_array/plots/correlation-rl.pdf
# data/shape-randst/0_25-18_60-28_60-94_00-15_00-full_array/plots/correlation-rl.pdf
# data/shape-randst/0_25-18_50-18_50-95_00-15_00-full_array/plots/correlation-rl.pdf
# data/shape-stad/0_25-2_92-1_18-0_00-15_00-full_array/plots/correlation-rl.pdf
# data/shape-stad/0_25-2_35-1_32-0_00-15_00-full_array/plots/correlation-rl.pdf

