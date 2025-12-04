#!/bin/bash 

#OAR -l host=16/cpu=2,walltime=6:00:00
#OAR -O results/mitm_OAR_%jobid%.out
#OAR -E results/mitm_OAR_%jobid%.err
#OAR -p paradoxe

N=30


GROUPS_COUNT_FILL=(16 8 6 4 2)
GROUPS_COUNT_PROBE=(16 8 6 4 2)
BUFFER_SIZE=(65536 32768 16384 8192 4096 2048 1024 512 256 128)
BSEND_AMOUNT=(1000)


gcc -Wall -o make_header make_header.c


for fill in "${GROUPS_COUNT_FILL[@]}"; do
    for probe in "${GROUPS_COUNT_PROBE[@]}"; do
        for buffer in "${BUFFER_SIZE[@]}"; do
            for bsend in "${BSEND_AMOUNT[@]}"; do
                echo "Running with fill_groups=$fill, probe_groups=$probe, buffer_size=$buffer, bsend_amount=$bsend"
                ./make_header --fill-groups $fill --probe-groups $probe --buffer-size $buffer --bsend-amount $bsend

                # We compile and run on the compute nodes because of the -march=native flag
                mpicc -O3 -march=native -Wall -o mitm_mpi_grouped mitm_mpi_grouped.c -fopenmp

                mpiexec \
                    --mca plm_rsh_agent oarsh \
                    --mca pml ob1 \
                    --mca btl ^openib \
                    --bind-to none \
                    --map-by ppr:1:node \
                    --hostfile $OAR_NODEFILE \
                    ./mitm_mpi_grouped --n $N --online
            done
        done
    done
done