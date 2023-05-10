#!/bin/bash

sys='sysA sysB evbless'

for j in ${sys[@]}; do
    #for i in fep_0*_${j}.log; do
        if grep "Finished mdrun on rank" fep_050_${j}.log >& /dev/null; then
            echo "$j:  - OK"
        else
	    echo "$j:  - FAILED"
        fi
    #done
done
