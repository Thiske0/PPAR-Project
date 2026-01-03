# Pascal Kessler, Mathis Poppe

# ================= GLOBAL SETTINGS =================

set terminal pngcairo size 900,600 enhanced font "Helvetica,12"
set logscale x
set logscale y
set grid
set key top left
set xlabel "Number of hosts"
set ylabel "Time (s)"

# Column mapping
COL_N     = 1
COL_FILL  = 2
COL_PROBE = 3
COL_TOTAL = 4

# Line styles for algorithms
set style line 1 lw 2 pt 7   # MPI
set style line 2 lw 2 pt 5   # MPI no comm
set style line 3 lw 2 pt 9   # NUMA
set style line 4 lw 2 pt 11  # NUMA no comm
set style line 5 lw 2 pt 13
set style line 6 lw 2 pt 15

set output "results/final results/performance_fill_multiple.png"
set title "Fill - Multiple host performance (N=30)"

plot \
  "results/final results/h 16 n 30 mpi reduced.txt" using COL_N:COL_FILL  with linespoints ls 1 title "MPI", \
  "results/final results/h 16 n 30 mpi no comm reduced.txt" using COL_N:COL_FILL  with linespoints ls 2 title "MPI no comm", \
  "results/final results/h 16 n 30 numa reduced.txt" using COL_N:COL_FILL   with linespoints ls 3 title "NUMA", \
  "results/final results/h 16 n 30 numa no comm reduced.txt" using COL_N:COL_FILL   with linespoints ls 4 title "NUMA no comm"

set output "results/final results/performance_probe_multiple.png"
set title "Probe - Multiple host performance (N=30)"

plot \
  "results/final results/h 16 n 30 mpi reduced.txt" using COL_N:COL_FILL  with linespoints ls 1 title "MPI", \
  "results/final results/h 16 n 30 mpi no comm reduced.txt" using COL_N:COL_PROBE  with linespoints ls 2 title "MPI no comm", \
  "results/final results/h 16 n 30 numa reduced.txt" using COL_N:COL_PROBE   with linespoints ls 3 title "NUMA", \
  "results/final results/h 16 n 30 numa no comm reduced.txt" using COL_N:COL_PROBE   with linespoints ls 4 title "NUMA no comm"

set output "results/final results/performance_total_multiple.png"
set title "Total - Multiple host performance (N=30)"

plot \
  "results/final results/h 16 n 30 mpi reduced.txt" using COL_N:COL_FILL  with linespoints ls 1 title "MPI", \
  "results/final results/h 16 n 30 mpi no comm reduced.txt" using COL_N:COL_TOTAL  with linespoints ls 2 title "MPI no comm", \
  "results/final results/h 16 n 30 numa reduced.txt" using COL_N:COL_TOTAL   with linespoints ls 3 title "NUMA", \
  "results/final results/h 16 n 30 numa no comm reduced.txt" using COL_N:COL_TOTAL   with linespoints ls 4 title "NUMA no comm"


set xlabel "Number of groups"

set output "results/final results/performance_communication_mpi.png"
set title "Communication performance (N=30, 16 hosts)"

plot \
  "results/final results/h 16 n 30 mpi large groups reduced.txt" using COL_N:COL_FILL  with linespoints ls 1 title "MPI Fill", \
  "results/final results/h 16 n 30 mpi large groups reduced.txt" using COL_N:COL_PROBE  with linespoints ls 2 title "MPI Probe", \
  "results/final results/h 16 n 30 mpi large groups reduced.txt" using COL_N:COL_TOTAL   with linespoints ls 3 title "MPI Total", \

set output "results/final results/performance_communication_numa.png"
set title "Communication performance (N=30, 16 hosts)"

plot \
  "results/final results/h 16 n 30 numa large groups reduced.txt" using COL_N:COL_FILL  with linespoints ls 1 title "NUMA Fill", \
  "results/final results/h 16 n 30 numa large groups reduced.txt" using COL_N:COL_PROBE  with linespoints ls 2 title "NUMA Probe", \
  "results/final results/h 16 n 30 numa large groups reduced.txt" using COL_N:COL_TOTAL   with linespoints ls 3 title "NUMA Total", \
