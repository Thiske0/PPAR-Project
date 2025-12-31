#!/bin/bash 

#OAR -l host=30/cpu=2,walltime=13:30:00
#OAR -O results/mitm_OAR_%jobid%.out
#OAR -E results/mitm_OAR_%jobid%.err
#OAR -p paradoxe

N=42
REDUCE=3

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

./make_header --fill-groups $GROUPS_COUNT_FILL --probe-groups $GROUPS_COUNT_PROBE --block-size $BLOCK_SIZE --buffer-size $BUFFER_SIZE --queue-size $QUEUE_SIZE --bsend-amount $BSEND_AMOUNT --prefill-buffer-size $PREFILL_BUFFER_SIZE

# We compile and run on the compute nodes because of the -march=native flag
mpicc -O3 -march=native -Wall -o mitm_numa mitm_numa_no_comm.c -fopenmp -lnuma

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

mpiexec \
    --mca plm_rsh_agent oarsh \
    --mca pml ob1 \
    --mca btl ^openib \
    --rankfile $tmp_rankfile \
    ./mitm_numa --n $N --online --reduce $REDUCE