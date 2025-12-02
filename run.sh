#!/bin/bash 

#OAR -l host=2,walltime=0:10:00
#OAR -O mandel_OAR_%jobid%.out
#OAR -E mandel_OAR_%jobid%.err

# We compile and run on the compute nodes because of the -march=native flag
mpicc -O3 -march=native -Wall -o mitm_mpi_grouped mitm_mpi_grouped.c -fopenmp
mpiexec -ppn 1 --hostfile $OAR_NODEFILE ./mitm_mpi_grouped --n 20 --online