#!/bin/bash


srun --job-name="single-im-cr p 3.00 0.50 full 500.00 4000.00" python pyplots/single-image-creation.py p 3.00 0.50 0.00 0.00 15.00 full_array 300.00 360.00 &
srun --job-name="single-im-cr p 3.00 0.50 exact 260.00 320.00" python pyplots/single-image-creation.py p 3.00 0.50 0.00 0.00 15.00 exact 260.00 320.00
