#!/bin/bash

sys='sysA sysB evbless'

for j in ${sys[@]}; do
    for i in fep_*_${j}.log; do
        if grep "Finished mdrun on rank" ${i} >& /dev/null; then
            echo "${i}:  - OK"
        else
	    echo "${i}:  - FAILED"
        fi
    done
    echo
done
