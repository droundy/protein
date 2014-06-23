#!/bin/bash

echo starting !
srun python density-correlation.py p 3.00 0.50 0.00 0.00 15.00 0.00 10.0
srun python density-correlation.py p 2.00 0.50 0.00 0.00 15.00 0.00 10.0
srun python density-correlation.py p 4.00 0.50 0.00 0.00 15.00 0.00 10.0
srun python density-correlation.py p 9.00 0.50 0.00 0.00 15.00 0.00 10.0
srun python density-correlation.py p 0.30 0.30 0.00 0.00 15.00 0.00 10.0
srun python density-correlation.py p 0.45 0.45 0.00 0.00 15.00 0.00 10.0

echo finished with the pills, now moving to the triangles
srun python density-correlation.py triangle 0.25 2.00 2.00 2.00 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 2.33 2.33 2.33 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 2.66 2.66 2.66 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 3.00 3.00 3.00 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 3.33 3.33 3.33 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 3.66 3.66 3.66 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 4.00 4.00 4.00 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 6.00 6.00 6.00 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 8.00 8.00 8.00 15.00 0.00 10.0
srun python density-correlation.py triangle 0.25 9.00 9.00 9.00 15.00 0.00 10.0
