#!/bin/bash
echo starting!

srun python pyplots/box_plot.py triangle 0.25 2.00 2.00 2.00 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 2.33 2.33 2.33 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 2.66 2.66 2.66 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 3.00 3.00 3.00 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 3.33 3.33 3.33 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 3.66 3.66 3.66 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 4.00 4.00 4.00 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 6.00 6.00 6.00 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 7.00 7.00 7.00 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 8.00 8.00 8.00 15.00 full_array 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 9.00 9.00 9.00 15.00 full_array 0.00 10.00


srun python pyplots/box_plot.py triangle 0.25 2.00 2.00 2.00 15.00 exact 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 2.33 2.33 2.33 15.00 exact 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 2.66 2.66 2.66 15.00 exact 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 3.00 3.00 3.00 15.00 exact 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 4.00 4.00 4.00 15.00 exact 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 6.00 6.00 6.00 15.00 exact 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 7.00 7.00 7.00 15.00 exact 0.00 10.00
srun python pyplots/box_plot.py triangle 0.25 8.00 8.00 8.00 15.00 exact 0.00 10.00
