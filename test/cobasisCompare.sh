#!/bin/bash

#tests difference of Basil and BasilP output on a set of examples 
# for each command line argument $ARG, there should be a file input/$ARG.txt of 
# Basil input, output/$ARG.bas.out of Basil output, and output/$ARG.bas.2.out 
# of BasilP output. This script will print the difference between the cobasis 
# lists of the two

for ARG in "$@"
do
	./basilCobasisRip.pl < "output/${ARG}.bas.2.out" > "isomorphTest.new.cob"
	./basilCobasisRip.pl < "output/${ARG}.bas.out" > "isomorphTest.old.cob"
	./ine2gapGrp.pl < "input/${ARG}.txt" > "isomorphTest-grp.gap"
	./gap.sh IsomorphicCompare.gap | ./catName.pl ":${ARG}:"
done;
