#!/bin/bash 

#OAR -l host=2/cpu=2,walltime=0:10:00
#OAR -O results/mitm_OAR_%jobid%.out
#OAR -E results/mitm_OAR_%jobid%.err
#OAR -p paradoxe

N=30


GROUPS_COUNT_FILL=1
GROUPS_COUNT_PROBE=1
BUFFER_SIZE=8192
BSEND_AMOUNT=1000

gcc -Wall -o make_header make_header.c

./make_header --fill-groups $GROUPS_COUNT_FILL --probe-groups $GROUPS_COUNT_PROBE --buffer-size $BUFFER_SIZE --bsend-amount $BSEND_AMOUNT

# We compile and run on the compute nodes because of the -march=native flag
mpicc -O3 -march=native -Wall -o mitm_mpi_grouped mitm_mpi_grouped.c -fopenmp

mpiexec \
    --mca pml ob1 \
    --mca btl ^openib \
    --bind-to none \
    --map-by ppr:1:node \
    --hostfile $OAR_NODEFILE \
    ./mitm_mpi_grouped --n $N --online