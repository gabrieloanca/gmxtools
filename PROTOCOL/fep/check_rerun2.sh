#!/bin/bash

sys='sysA sysB evbless'

for j in ${sys[@]}; do
    for i in fep_0*_${j}.log; do
        if grep "Finished mdrun on rank" ${i}.log >& /dev/null; then
            echo "$i:  - OK"
        else
	    echo "$i:  - FAILED"
        fi
    done
done
