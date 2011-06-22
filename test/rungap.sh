#!/bin/bash

if [ -z $1 ]; then
	gap -b -q -e 
else
	gap -b -q -e < $1
fi;
