#!/bin/bash

file=$1
#sed "s/nsttcouple               = -1        ; -1 = nstlist/nsttcouple               = -1        ; -1 = 10 or less (default)/" < ${file} > ${file}.tmp
#sed "s/nstpcouple               = -1        ; -1 = nstlist/nstpcouple               = -1        ; -1 = 10 or less (default)/" < ${file} > ${file}.tmp
sed "s/pcoupl                   = C-rescale/pcoupl                   = no        ; C-rescale/" < ${file} > ${file}.tmp



#sed "s///" < ${file} > ${file}.tmp
#sed -i 's///' ${file}
#for j in *.tmp; do mv $j `basename $j .tmp`; done
#for i in *.mdp; do ./fixfile.sh $i; done


