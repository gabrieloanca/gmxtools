#!/bin/bash

par1=$1

sys='sysA sysB evbless'

for k in ${sys[@]}; do
    for l in {000..050}; do
        file=fep_${l}_${k};
        if grep "Finished mdrun on rank" ${file}.log >& /dev/null; then
            if find . -name ${file}.edr >& /dev/null; then
                echo "${par1}/${file}:  - OK"
            else
                echo "${par1}/${file}:  - FAILED"
            fi
        else
	    echo "${par1}/${file}:  - FAILED"
        fi
    done
    echo
done
