#!/bin/bash

# srun --job-name="sim p 3.00 0.50 15.00 full" ./protein_microscopy p 3.00 0.50 0.00 0.00 15.00 1.00 full_array -dump &
# srun --job-name="sim p 3.00 0.50 15.00 exact" ./protein_microscopy p 3.00 0.50 0.00 0.00 15.00 1.00 exact -dump &

# srun --job-name="sim randst 0.25 8.00 8.00 95.00 15.00 0.5 exact" ./protein_microscopy randst 0.25 8.00 8.00 95.00 15.00 0.5 exact -dump &
# srun --job-name="sim randst 0.25 10.00 17.00 94.00 15.00 0.5 exact" ./protein_microscopy randst 0.25 10.00 17.00 94.00 15.00 0.5 exact -dump &

# srun --job-name="sim randst 0.25 18.50 18.50 95.00 15.00 1.35 full" ./protein_microscopy randst 0.25 18.50 18.50 95.00 15.00 1.35 full_array -dump &
# srun --job-name="sim randst 0.25 18.50 18.50 95.00 15.00 1.35 exact" ./protein_microscopy randst 0.25 18.50 18.50 95.00 15.00 1.35 exact -dump &

# srun --job-name="sim randst 0.25 18.60 28.60 94.00 15.00 1.53 full" ./protein_microscopy randst 0.25 18.60 28.60 94.00 15.00 1.53 full_array -dump &
# srun --job-name="sim randst 0.25 18.60 28.60 94.00 15.00 1.53 exact" ./protein_microscopy randst 0.25 18.60 28.60 94.00 15.00 1.53 exact -dump &

# srun --job-name="sim stad 0.25 5.00 1.00 15.00 full" ./protein_microscopy stad 0.25 5.00 1.00 0.00 15.00 1.0 full_array -dump &
# srun --job-name="sim stad 0.25 5.00 1.00 15.00 exact" ./protein_microscopy stad 0.25 5.00 1.00 0.00 15.00 1.0 exact -dump

# srun --job-name="sim stad 0.25 3.00 0.50 15.00 full" ./protein_microscopy stad 0.25 3.00 0.50 0.00 15.00 1.0 full_array -dump &
# srun --job-name="sim stad 0.25 3.00 0.50 15.00 exact" ./protein_microscopy stad 0.25 3.00 0.50 0.00 15.00 1.0 exact -dump &

# srun --job-name="sim stad 0.25 3.00 1.18 15.00 full" ./protein_microscopy stad 0.25 3.00 1.18 0.00 15.00 1.0 full_array -dump &
# srun --job-name="sim stad 0.25 3.00 1.18 15.00 exact" ./protein_microscopy stad 0.25 3.00 1.18 0.00 15.00 1.0 exact -dump &

# srun --job-name="sim stad 0.25 2.35 1.32 15.00 full" ./protein_microscopy stad 0.25 2.35 1.32 0.00 15.00 1.0 full_array -dump &
# srun --job-name="sim stad 0.25 2.35 1.32 15.00 exact" ./protein_microscopy stad 0.25 2.35 1.32 0.00 15.00 1.0 exact -dump &

# srun --job-name="sim stad 0.25 2.92 1.18 15.00 full" ./protein_microscopy stad 0.25 2.92 1.18 0.00 15.00 1.0 full_array -dump &
# srun --job-name="sim stad 0.25 2.92 1.18 15.00 exact" ./protein_microscopy stad 0.25 2.92 1.18 0.00 15.00 1.0 exact -dump &

##################################

# srun --mem-per-cpu=7000 --job-name="hi sim stad 0.25 2.91 1.18 15.00 full" ./protein_microscopy stad 0.25 2.91 1.18 0.00 15.00 1.0 full_array -dump -hires &
# srun --mem-per-cpu=7000 --job-name="hi sim stad 0.25 2.91 1.18 15.00 exact" ./protein_microscopy stad 0.25 2.91 1.18 0.00 15.00 1.0 exact -dump -hires &

# srun --mem-per-cpu=7000 --job-name="hi sim stad 0.25 2.34 1.32 15.00 full" ./protein_microscopy stad 0.25 2.34 1.32 0.00 15.00 1.0 full_array -dump -hires  &
# srun --mem-per-cpu=7000 --job-name="hi sim stad 0.25 2.34 1.32 15.00 exact" ./protein_microscopy stad 0.25 2.34 1.32 0.00 15.00 1.0 exact -dump -hires &

# srun --mem-per-cpu=7000 --job-name="hi sim randst 0.25 18.51 18.50 95.00 15.00 1.35 full" ./protein_microscopy randst 0.25 18.51 18.50 95.00 15.00 1.35 full_array -dump -hires &
# srun --mem-per-cpu=7000 --job-name="hi sim randst 0.25 18.51 18.50 95.00 15.00 1.35 exact" ./protein_microscopy randst 0.25 18.51 18.50 95.00 15.00 1.35 exact -dump -hires &
# srun --mem-per-cpu=7000 --job-name="hi sim randst 0.25 18.61 28.60 94.00 15.00 1.53 full" ./protein_microscopy randst 0.25 18.61 28.60 94.00 15.00 1.53 full_array -dump -hires &
# srun --mem-per-cpu=7000 --job-name="hi sim randst 0.25 18.61 28.60 94.00 15.00 1.53 exact" ./protein_microscopy randst 0.25 18.61 28.60 94.00 15.00 1.53 exact -dump -hires &

# srun --job-name="sim p 3.00 0.50 15.00 hi full" ./protein_microscopy p 3.00 0.50 0.00 0.00 15.00 1.00 full_array -dump -hires &
# srun --job-name="sim p 3.00 0.50 15.00 hi exact" ./protein_microscopy p 3.00 0.50 0.00 0.00 15.00 1.00 exact -dump -hires &
