#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --time=00:15:00
#SBATCH --output=/WeakScalingFDSTAG_2cores.txt

# this is to perform an isoviscous weak scaling test on rosa, with 48^3 per node

aprun -n 2 ../../bin/LaMEM -options_file options_isoviscous_rosa.txt
