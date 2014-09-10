#!/bin/bash

echo starting!
# srun python pyplots/movie-image-creation.py p 1.00 0.50 0.00 0.00 15.00 full_array 300.00 1000.00 &
# srun python pyplots/movie-image-creation.py p 2.00 0.50 0.00 0.00 15.00 full_array 300.00 1000.00 &

srun --job-name="mov-im-create p 3.00 0.50 full 500 1000" python pyplots/movie-image-creation.py p 3.00 0.50 0.00 0.00 15.00 full_array 500.00 1000.00 &
srun --job-name="mov-im-create p 3.00 0.50 exact 500 1000" python pyplots/movie-image-creation.py p 3.00 0.50 0.00 0.00 15.00 exact 500.00 1000.00

# srun python pyplots/movie-image-creation.py p 4.00 0.50 0.00 0.00 15.00 full_array 300.00 1000.00 &

# srun python pyplots/movie-image-creation.py stad 0.25 6.00 2.00 0.00 15.00 full_array 300.00 1000.00 &

# srun python pyplots/movie-image-creation.py stad 0.25 3.00 0.50 0.00 15.00 full_array 500.00 1000.00 &
# srun python pyplots/movie-image-creation.py stad 0.25 3.00 0.50 0.00 15.00 exact 500.00 1000.00 &

# srun python pyplots/movie-image-creation.py stad 0.25 3.00 1.18 0.00 15.00 full_array 500.00 1000.00 &
# #srun python pyplots/movie-image-creation.py stad 0.25 3.00 1.18 0.00 15.00 exact 500.00 1000.00 &

# srun python pyplots/movie-image-creation.py randst 0.25 18.50 18.50 95.00 15.00 full_array 500.00 1000.00 &
# srun python pyplots/movie-image-creation.py randst 0.25 18.50 18.50 95.00 15.00 exact 500.00 1000.00 &

# srun --job-name="mov-im-create randst 18.60 28.60 94.00 full 500 1100" python pyplots/movie-image-creation.py randst 0.25 18.60 28.60 94.00 15.00 full_array 500.00 1000.00 &
# srun --job-name="mov-im-create randst 18.60 28.60 94.00 exact 500 1100" python pyplots/movie-image-creation.py randst 0.25 18.60 28.60 94.00 15.00 exact 500.00 1500.00 &

# srun --job-name="mov-make randst 18.60 28.60 94.00 full 500 1100" python pyplots/movie-maker.py randst 0.25 18.60 28.60 94.00 15.00 full_array 500.00 1000.00 &
# srun --job-name="mov-make randst 18.60 28.60 94.00 exact 500 1100" python pyplots/movie-maker.py randst 0.25 18.60 28.60 94.00 15.00 exact 500.00 1500.00 &

# srun --job-name="mov-im-create stad 2.35 1.32 full 500.00 1500.00" python pyplots/movie-image-creation.py stad 0.25 2.35 1.32 0.00 15.00 full_array 500.00 1100.00 &
# srun --job-name="mov-im-create stad 2.35 1.32 exact 500.00 800.00" python pyplots/movie-image-creation.py stad 0.25 2.35 1.32 0.00 15.00 exact 500.00 810.00 &

# srun --job-name="mov-make stad 2.35 1.32 full 500.00 1500.00" python pyplots/movie-maker.py stad 0.25 2.35 1.32 0.00 15.00 full_array 500.00 1100.00 &
# srun --job-name="mov-make stad 2.35 1.32 exact 500.00 800.00" python pyplots/movie-maker.py stad 0.25 2.35 1.32 0.00 15.00 exact 500.00 810.00 &

# srun --job-name="mov-im-create stad 2.92 1.18 full 500.00 1500.00" python pyplots/movie-image-creation.py stad 0.25 2.92 1.18 0.00 15.00 full_array 500.00 1100.00 &
# srun --job-name="mov-im-create stad 2.92 1.18 exact 500.00 800.00" python pyplots/movie-image-creation.py stad 0.25 2.92 1.18 0.00 15.00 exact 500.00 810.00 &

# srun --job-name="mov-make stad 2.92 1.18 full 500.00 1500.00" python pyplots/movie-maker.py stad 0.25 2.92 1.18 0.00 15.00 full_array 500.00 1100.00 &
# srun --job-name="mov-make stad 2.92 1.18 exact 500.00 800.00" python pyplots/movie-maker.py stad 0.25 2.92 1.18 0.00 15.00 exact 500.00 810.00 &


# srun python pyplots/movie-image-creation.py stad 0.25 3.00 0.50 0.00 15.00 full_array 500.00 1000.00 &
# srun python pyplots/movie-image-creation.py stad 0.25 3.00 0.50 0.00 15.00 exact 500.00 1000.00 &


# arg_set = ["randst/0.25-18.50-18.50-95.00-15.00-exact",
#            "randst/0.25-18.50-18.50-95.00-15.00-full_array",
#            "randst/0.25-18.60-28.60-94.00-15.00-exact",
#            "randst/0.25-18.60-28.60-94.00-15.00-full_array",
#            "stad/0.25-2.35-1.32-0.00-15.00-full_array",
#            "stad/0.25-2.35-1.32-0.00-15.00-exact",
#            "stad/0.25-2.92-1.18-0.00-15.00-full_array",
#            "stad/0.25-2.92-1.18-0.00-15.00-exact"]

#bound_times = [500,1000,500,980,500,1000,500,1000]
