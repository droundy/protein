#!/bin/bash

echo starting!
srun python pyplots/density_movie.py p 0.45 0.45 0.00 0.00 15.00 full_array 0.00 900.00 &
srun python pyplots/density_movie.py p 2.00 0.50 0.00 0.00 15.00 full_array 0.00 900.00 &
srun python pyplots/density_movie.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00 900.00 &
srun python pyplots/density_movie.py p 4.00 0.50 0.00 0.00 15.00 full_array 0.00 900.00 &
srun python pyplots/density_movie.py p 9.00 0.50 0.00 0.00 15.00 full_array 0.00 900.00 &
