#!/bin/bash

#=========================================================================================
# Script to lunch parallel MPI application on Linux under valgrind with MPI wrappers
#
# All errors and screen output will be redirected to the file: 
#      outfile_name.out
#
# Valgrind errors & warnings will be written to a separate file for each MPI task: 
#      outfile_name.pid.xml (pid - system process ID of MPI task, assigned automatically)
# 
# To visualize memory check summary in human readable format use Valkyrie GUI.
#
# Invocation pattern:
#
# ./ValgrindCheck.sh num_proc "exec_path prog_args" outfile_name
#
#=========================================================================================

MPIWRAP_DEBUG=quiet \
LD_PRELOAD=/opt/valgrind/lib/valgrind/libmpiwrap-amd64-linux.so \
mpirun -np $1 \
valgrind -v \
--leak-check=full \
--track-origins=yes \
--show-reachable=yes \
--xml=yes \
--xml-file=$3_%p.xml \
--child-silent-after-fork=yes \
-q \
$2 > $3.out 2>&1

exit
