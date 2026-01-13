#!/bin/bash

mpgl=$1
K=$2
r=$3

echo
echo "input file:" $mpgl
echo "K:" $K
echo "replication:" $r
echo
echo "-------------entropy start-------------"

entropy -i $mpgl -m 1 -n ploidy_inds.txt \
        -k $K -q qk${K}inds.txt \
		-l 10000 -b 1000 -t 5 -D 1\
		-o Result/mcmcoutchain${r}.k${K}.hdf5