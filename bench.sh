#!/bin/bash 

#OAR -l host=16/cpu=2,walltime=6:00:00
#OAR -O results/mitm_OAR_%jobid%.out
#OAR -E results/mitm_OAR_%jobid%.err
#OAR -p paradoxe

N=30
REDUCE=0
PROG=mitm_numa_no_comm

THREADS=(104)
HOSTS=(16)

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

gcc -Wall -o make_header make_header.c

# create temporary file
tmp_rankfile=$(mktemp)

# Get unique nodes, sorted numerically
mapfile -t NODES < <(
    uniq "$OAR_NODEFILE" | sort -V
)

rank=0
> "$tmp_rankfile"

for node in "${NODES[@]}"; do
    echo "rank $rank=$node slot=0-51" >> "$tmp_rankfile" # use all the slots on the node
    rank=$((rank + 1))
done

# We compile and run on the compute nodes because of the -march=native flag

for host in "${HOSTS[@]}"; do
    reduced_rank_file=$(mktemp)
    head -n $host $tmp_rankfile > $reduced_rank_file

    ./make_header --fill-groups $GROUPS_COUNT_FILL --probe-groups $GROUPS_COUNT_PROBE --block-size $BLOCK_SIZE --buffer-size $BUFFER_SIZE --queue-size $QUEUE_SIZE --bsend-amount $BSEND_AMOUNT --prefill-buffer-size $PREFILL_BUFFER_SIZE
    mpicc -O3 -march=native -Wall -o $PROG $PROG.c -fopenmp -lnuma

    for threads in "${THREADS[@]}"; do

        export OMP_NUM_THREADS=$threads
        mpiexec \
            --mca plm_rsh_agent oarsh \
            --mca pml ob1 \
            --mca btl ^openib \
            --rankfile $reduced_rank_file \
            ./$PROG --n $N --online
    done
done