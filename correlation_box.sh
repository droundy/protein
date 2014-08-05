#!/bin/bash

echo starting !
# srun python density-correlation.py p 3.00 0.50 0.00 0.00 15.00 0.00 10.0 &
# srun python density-correlation.py p 2.00 0.50 0.00 0.00 15.00 0.00 10.0 &
# srun python density-correlation.py p 4.00 0.50 0.00 0.00 15.00 0.00 10.0 &
# srun python density-correlation.py p 9.00 0.50 0.00 0.00 15.00 0.00 10.0 &
# srun python density-correlation.py p 0.30 0.30 0.00 0.00 15.00 0.00 10.0 &
# srun python density-correlation.py p 0.45 0.45 0.00 0.00 15.00 0.00 10.0 &

# echo finished with the pills, now moving to the triangles
# srun python density-correlation.py triangle 0.25 2.00 2.00 2.00 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 2.33 2.33 2.33 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 2.66 2.66 2.66 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 3.00 3.00 3.00 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 3.33 3.33 3.33 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 3.66 3.66 3.66 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 4.00 4.00 4.00 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 6.00 6.00 6.00 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 8.00 8.00 8.00 15.00 0.00 10.0 &
# srun python density-correlation.py triangle 0.25 9.00 9.00 9.00 15.00 0.00 10.0 &



srun python density-correlation.py randst 0.25 10.00 17.00 94.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 12.00 19.00 94.00 15.00 0.00 10.00 &

srun python density-correlation.py randst 0.25 8.00 8.00 95.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 10.00 10.00 95.00 15.00 0.00 10.00 &

srun python density-correlation.py randst 0.25 15.00 15.00 95.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 17.00 17.00 95.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 18.50 18.50 95.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 21.00 21.00 95.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 25.00 25.00 95.00 15.00 0.00 10.00 &

srun python density-correlation.py randst 0.25 14.00 24.00 94.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 16.00 26.00 94.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 18.50 28.00 94.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 18.60 28.60 94.00 15.00 0.00 10.00 &
srun python density-correlation.py randst 0.25 23.00 32.00 94.00 15.00 0.00 10.00 &

srun python density-correlation.py stad 0.25 5.00 2.50 0.00 15.00 0.00 10.00 &
srun python density-correlation.py stad 0.25 6.00 2.00 0.00 15.00 0.00 10.00 &

srun python density-correlation.py p 4.00 0.50 0.00 0.00 15.00 0.00 10.00 &
