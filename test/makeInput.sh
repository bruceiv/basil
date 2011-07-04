#!/bin/bash

#generates Basil and Symbal input files.
# PRE-CONDITIONS:
# 1. The constraint matrix is stored in GAP format in nautyMat.gap, named M
# 2. This script must be called with one argument, the name of the test case

if [ -z $1 ]; then
	echo -e "usage:\n  makeinput.sh <input_name>"
else
	# use Nauty to generate the automorphism group for the matrix
	./gap.sh nautyGen.gap > ${1}"-grp.gap"
	# convert the generated group to Basil format
	./gapGrp2ine.pl < ${1}"-grp.gap" > ${1}"-grp.txt"
	# convert the input matrix to Basil format
	./gap2ine.pl < nautyMat.gap > ${1}".ine"
	# combine the matrix and group into a Basil input file
	echo -e "${1}\nH-representation" \
		| cat - ${1}".ine" ${1}"-grp.txt" \
		> ${1}".txt"
	# convert this Basil input file to a Symbal file
	./ine2gap.pl < ${1}".txt" > ${1}".gap"
	# remove intermediate files
	rm ${1}"-grp.gap" ${1}"-grp.txt" ${1}".ine"
fi;