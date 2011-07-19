#!/bin/bash

#tests difference of Basil and Symbal output on a set of examples
# for each command line argument $ARG, there should be a file $ARG.bas.out of 
# Basil output, and $ARG.sym.out of Symbal output. This script will print the 
# difference between the cobasis lists of the two

for ARG in "$@"
do
	./basilCobasisRip.pl < "${ARG}.bas.out" > "isomorphTest.bas.cob"
	./symbalCobasisRip.pl < "${ARG}.sym.out" > "isomorphTest.sym.cob"
	cp "${ARG}-grp.gap" "isomorphTest-grp.gap"
	./gap.sh IsomorphicDiff.gap | ./catName.pl ":${ARG}:"
done;
