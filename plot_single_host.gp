# Pascal Kessler, Mathis Poppe

# ================= GLOBAL SETTINGS =================

set terminal pngcairo size 900,600 enhanced font "Helvetica,12"
set logscale x
set logscale y
set grid
set key top left
set xlabel "Number of threads"
set ylabel "Time (s)"

# Column mapping
COL_N     = 1
COL_FILL  = 2
COL_PROBE = 3
COL_TOTAL = 4

# Line styles for algorithms
set style line 1 lw 2 pt 7   # Sequential
set style line 2 lw 2 pt 5   # Vectorized
set style line 3 lw 2 pt 9   # OpenMP
set style line 4 lw 2 pt 11  # NUMA

set output "results/final results/performance_fill_single.png"
set title "Fill - Single host performance (N=27)"

plot \
  "results/final results/seq n 27 reduced.txt" using COL_N:COL_FILL  with linespoints ls 1 title "Sequential", \
  "results/final results/vectorized n 27 reduced.txt" using COL_N:COL_FILL  with linespoints ls 2 title "Vectorized", \
  "results/final results/omp n 27 reduced.txt" using COL_N:COL_FILL  with linespoints ls 3 title "OpenMP", \
  "results/final results/numa h 1 n 27 reduced.txt" using COL_N:COL_FILL   with linespoints ls 4 title "NUMA"

set output "results/final results/performance_probe_single.png"
set title "Probe - Single host performance (N=27)"

plot \
  "results/final results/seq n 27 reduced.txt" using COL_N:COL_PROBE  with linespoints ls 1 title "Sequential", \
  "results/final results/vectorized n 27 reduced.txt" using COL_N:COL_PROBE  with linespoints ls 2 title "Vectorized", \
  "results/final results/omp n 27 reduced.txt" using COL_N:COL_PROBE  with linespoints ls 3 title "OpenMP", \
  "results/final results/numa h 1 n 27 reduced.txt" using COL_N:COL_PROBE   with linespoints ls 4 title "NUMA"

set output "results/final results/performance_total_single.png"
set title "Total - Single host performance (N=27)"

plot \
  "results/final results/seq n 27 reduced.txt" using COL_N:COL_TOTAL  with linespoints ls 1 title "Sequential", \
  "results/final results/vectorized n 27 reduced.txt" using COL_N:COL_TOTAL  with linespoints ls 2 title "Vectorized", \
  "results/final results/omp n 27 reduced.txt" using COL_N:COL_TOTAL  with linespoints ls 3 title "OpenMP", \
  "results/final results/numa h 1 n 27 reduced.txt" using COL_N:COL_TOTAL   with linespoints ls 4 title "NUMA"