#!/bin/bash

srun --job-name="sim p 3.00 0.50 15.00 full" ./protein_microscopy p 3.00 0.50 0.00 0.00 15.00 1.00 full_array -dump &
srun --job-name="sim p 3.00 0.50 15.00 exact" ./protein_microscopy p 3.00 0.50 0.00 0.00 15.00 1.00 exact -dump &

srun --job-name="sim randst 0.25 18.50 18.50 95.00 15.00 1.35 full" ./protein_microscopy randst 0.25 18.50 18.50 95.00 15.00 1.35 full_array -dump &
srun --job-name="sim randst 0.25 18.50 18.50 95.00 15.00 1.35 exact" ./protein_microscopy randst 0.25 18.50 18.50 95.00 15.00 1.35 exact -dump &

srun --job-name="sim randst 0.25 18.60 28.60 94.00 15.00 1.53 full" ./protein_microscopy randst 0.25 18.60 28.60 94.00 15.00 1.53 full_array -dump &
srun --job-name="sim randst 0.25 18.60 28.60 94.00 15.00 1.53 exact" ./protein_microscopy randst 0.25 18.60 28.60 94.00 15.00 1.53 exact -dump &

srun --job-name="sim stad 0.25 5.50 1.00 15.00 full" ./protein_microscopy stad 0.25 5.50 1.00 0.00 15.00 1.0 full_array -dump &
srun --job-name="sim stad 0.25 5.50 1.00 15.00 exact" ./protein_microscopy stad 0.25 5.50 1.00 0.00 15.00 1.0 exact -dump &

srun --job-name="sim p 3.50 0.50 15.00 full" ./protein_microscopy p 3.50 0.50 0.00 0.00 15.00 1.0 full_array -dump &
srun --job-name="sim p 3.50 0.50 15.00 exact" ./protein_microscopy p 3.50 0.50 0.00 0.00 15.00 1.0 exact -dump &
