#!/bin/bash

# Create archives for all timesteps 
#
# --> Use it as follows
#
# 1.) ./archiveTS.sh
# 2.) ./archiveTS.sh -deldir
# 3.) ./archiveTS.sh -untar
#
# two possible options: 
# -deldir # removes directory after creating an archive
# -untar  # extract Timestep-archives
#
#
# Tobias Baumann, Mainz 2015


if [ "-untar" = "$1" ]
then
	for k in Timestep_*_*[0-9].tar
	do
		echo "> extract ${k}"
		tar -xf ${k}
	done
else
	for k in Timestep_*_*[0-9]
	do
		echo "${k}"
		tar -cf "${k}.tar" ${k}
		if [ "-deldir" = "$1" ]
		then
			echo "> delete ${k}"
			rm -rf ${k}
		fi
	done
fi
