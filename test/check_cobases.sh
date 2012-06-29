#/bin/sh

if [ -z $1 ]; then
	echo -e "usage:\n  check_cobases.sh <input_name>"
else
	./basil2oct.pl < $1 > bas_out.m
	octave -q cobase_check.m
	rm bas_out.m
fi;
