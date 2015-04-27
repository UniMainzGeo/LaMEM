#!/bin/bash

# Create archives for all timesteps 
# ./archiveTS.sh
#
# Tobias Baumann, Mainz 2015


for k in Timestep_*_*[0-9]
do
	tar -cf "${k}.tar" ${k}
done
