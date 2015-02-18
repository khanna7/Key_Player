#!/bin/bash

#SBATCH --ntasks=8

export R_LIBS=~/R_libs

module load R/3.1+intel-15.0 openmpi/1.8+intel-15.0
mpirun -np 1 Rscript bet-test-mpi.R
