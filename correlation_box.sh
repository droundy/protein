#!/bin/bash

echo starting !

# srun --job-name="corr randst 18.50 18.50 95.00 full_array 0.00 4500.00" python pyplots/density-correlation.py randst 0.25 18.50 18.50 95.00 15.00 full_array 0.00 4500.00 &
# srun --job-name="corr randst 18.50 18.50 95.00 exact 0.00 5000.00" python pyplots/density-correlation.py randst 0.25 18.50 18.50 95.00 15.00 exact 0.00 5000.00 &

# srun --job-name="corr randst 18.60 28.60 94.00 full_array 0.00 1020.00" python pyplots/density-correlation.py randst 0.25 18.60 28.60 94.00 15.00 full_array 0.00 1020.00 &
# srun --job-name="corr randst 18.60 28.60 94.00 exact 0.00 1800.00" python pyplots/density-correlation.py randst 0.25 18.60 28.60 94.00 15.00 exact 0.00 1800.00 &

srun --job-name="corr stad 2.35 1.32 full 0.00 10000.00" python pyplots/density-correlation.py stad 0.25 2.35 1.32 0.00 15.00 full_array 0.00 13500.00 &
srun --job-name="corr stad 2.35 1.32 exact 0.00 6000.00" python pyplots/density-correlation.py stad 0.25 2.35 1.32 0.00 15.00 exact 0.00 8400.00 &

# srun --job-name="corr stad 2.92 1.18 full 0.00 1410.00" python pyplots/density-correlation.py stad 0.25 2.92 1.18 0.00 15.00 full_array 0.00 1410.00 &
# srun --job-name="corr stad 2.92 1.18 exact 0.00 910.00" python pyplots/density-correlation.py stad 0.25 2.92 1.18 0.00 15.00 exact 0.00 910.00 &

# srun --job-name="corr p 3.00 0.50 full 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 6000.00 &
# srun --job-name="corr p 3.00 0.50 exact 0.00 6000.00" python pyplots/density-correlation.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00 6000.00

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
