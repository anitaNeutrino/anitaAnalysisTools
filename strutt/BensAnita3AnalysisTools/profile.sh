#!/bin/sh
##################################################################################################################
# Author: Ben Strutt
# Email: b.strutt.12@ucl.ac.uk
# Description: Systematically store the output of gperftools profiling
##################################################################################################################

# Developed on OS X, will require modification for Linux systems
# Can't get libtcmalloc to work on current system... so will only do CPUPROFILE for now

#set -x # echo commands as they are executed for debugging

outputDir="gperftoolsProfiles"
timeStamp=$(date "+%Y-%m-%d_%H-%M-%S" | sed -e 's/ /_/g')
prog=$1;

if [ -x $prog ]; then
    if [ ! -d $prog ]; then
	# Generates output directory if it doesn't exist
	if [ ! -d $outputDir ]; then
	    mkdir $outputDir
	fi;

	CPUPROFILE=/tmp/$prog"_"$timeStamp.prof ./$prog;
	pprof --pdf $prog /tmp/$prog"_"$timeStamp.prof > $outputDir/$prog"_"$timeStamp"_cpuProf.pdf"
    fi;
fi;
