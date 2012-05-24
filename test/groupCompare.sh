#!/bin/bash

#tests difference of Euclidean and Q-matrix symmetry groups
# for each command line argument $ARG, there should be a file input/$ARG.EG.txt 
# of Basil input with a Euclidean-generated group and input/$ARG.QG.txt of 
# Basil input with a Q-matrix-generated group. This script will print the 
# difference between the symmetry groups of the two

for ARG in "$@"
do
	./ine2gapGrp.pl < "input/${ARG}.EG.txt" > "e-grp.gap"
	./ine2gapGrp.pl < "input/${ARG}.QG.txt" > "q-grp.gap"
	./gap.sh GroupCompare.gap | ./catName.pl ":${ARG}:"
done;
rm e-grp.gap q-grp.gap
