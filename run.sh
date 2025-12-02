#!/bin/bash 

#OAR -l host=2/core=18,walltime=0:10:00
#OAR -O mitm_OAR_%jobid%.out
#OAR -E mitm_OAR_%jobid%.err

cat $OAR_NODEFILE

# We compile and run on the compute nodes because of the -march=native flag
mpicc -O3 -march=native -Wall -o mitm_mpi_grouped mitm_mpi_grouped.c -fopenmp
mpiexec  --map-by ppr:1:node --hostfile $OAR_NODEFILE ./mitm_mpi_grouped --n 20 --online