#!/bin/bash

# runs GAP without banner or prompts, in interpreter mode
# input file will be first argument if given, standard input otherwise
# output file will be second argument if given, standard output otherwise

if [ -z $1 ]; then
	gap -b -q -e 
else
	if [ -z $2 ]; then
		gap -b -q -e < $1
	else
		gap -b -q -e < $1 > $2
	fi;
fi;
