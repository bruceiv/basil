#!/bin/bash

#gets the order of a set of symmetry groups
# for each command line argument $ARG, there should be a file input/$ARG.txt 
# of Basil input. This script will print the order of the input group of that 
# input file.

for ARG in "$@"
do
	./ine2gapGrp.pl < "input/${ARG}.txt" > "bas-grp.gap"
	./gap.sh GroupOrder.gap | ./catName.pl ":${ARG}:"
done;
rm bas-grp.gap
