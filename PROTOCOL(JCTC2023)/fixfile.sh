#!/bin/bash

# Use it as:
# for i in *.mdp; do ./fixfile.sh $i; done

file=$1

sed 's/nsteps                   = 20000/nsteps                   = 10000/' < ${file} > ${file}.tmp 

#sed "s///" < ${file} > ${file}.tmp
#sed -i 's///' ${file}
#for j in *.tmp; do mv $j `basename $j .tmp`; done

