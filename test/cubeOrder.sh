#!/bin/bash

for i in {1..100}; do
	../basil -p --no-gram "cube/6-6-3-${i}.txt" | ./ine2gapGrp.pl > "bas-grp.gap"
	echo "6-6-3-${i}:" $(./gap.sh GroupOrder.gap)
done
for i in {1..100}; do
	../basil -p --no-gram "cube/7-6-3-${i}.txt" | ./ine2gapGrp.pl > "bas-grp.gap"
	echo "7-6-3-${i}:" $(./gap.sh GroupOrder.gap)
done
rm bas-grp.gap
