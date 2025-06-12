#!/bin/bash
set -e
mkdir -p build
g++ -O3 -fopenmp main.cpp physics.cpp -o riemann1d
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
./riemann1d

python3 plot_density.py
python3 visualize_mhd.py
