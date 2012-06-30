#!/bin/bash

#tests the validity of the isomorphisms for a set of Basil instances.
# for each command line argument $ARG, there should be a file input/$ARG.txt of 
# Basil input and output/$ARG.bas.out of Basil output. This script will print 
# each generator, and whether it passes the isomorphism test.

for ARG in "$@"
do
	./basil2oct.pl < "output/${ARG}.bas.out" > "img_mat.m"
	./basilGenRip.pl < "output/${ARG}.bas.out" > "img-gens.gap"
	./gap.sh CobImages.gap > "cob_imgs.m"
	octave -q isomorph_check.m | ./catName.pl ":${ARG}:"
done;
