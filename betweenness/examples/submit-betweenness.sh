#!/bin/bash

#SBATCH --ntasks=8

source /etc/profile

module load R/3.1+intel-15.0 openmpi/1.8+intel-15.0
mpirun -np 1 Rscript process-facebook-data.R
