#!/bin/bash 

#OAR -l host=16/cpu=2,walltime=6:00:00
#OAR -O results/mitm_OAR_%jobid%.out
#OAR -E results/mitm_OAR_%jobid%.err
#OAR -p paradoxe

N=27
PROG=mitm_numa_no_comm

# Optimal parameters for Paradoxe cluster for mitm_numa_no_comm
BLOCK_SIZE=131072

# Optimal parameters for Paradoxe cluster for mitm_numa
#BLOCK_SIZE=1024
QUEUE_SIZE=8192
BUFFER_SIZE=131072
PREFILL_BUFFER_SIZE=128

GROUPS_COUNT_FILL=2
GROUPS_COUNT_PROBE=2
BSEND_AMOUNT=20

set -euo pipefail
IFS=$'\n\t'

./make_header --fill-groups $GROUPS_COUNT_FILL --probe-groups $GROUPS_COUNT_PROBE --block-size $BLOCK_SIZE --buffer-size $BUFFER_SIZE --queue-size $QUEUE_SIZE --bsend-amount $BSEND_AMOUNT --prefill-buffer-size $PREFILL_BUFFER_SIZE


THREADS=(2 4 8 16 32 64 104)

gcc -Wall -o make_header make_header.c


for threads in "${THREADS[@]}"; do

    # We compile and run on the compute nodes because of the -march=native flag
    mpicc -O3 -march=native -Wall -o $PROG $PROG.c -fopenmp -lnuma

    export OMP_NUM_THREADS=$threads
    mpiexec \
        --mca plm_rsh_agent oarsh \
        --mca pml ob1 \
        --mca btl ^openib \
        --bind-to none \
        --map-by ppr:1:node \
        --hostfile $OAR_NODEFILE \
        ./$PROG --n $N --online
done