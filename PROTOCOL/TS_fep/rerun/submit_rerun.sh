#!/bin/bash

for k in {000..050}; do sbatch rerun.sh sysA fep_${k}; sleep 1; done
for k in {000..050}; do sbatch rerun.sh sysB fep_${k}; sleep 1; done
for k in {000..050}; do sbatch rerun.sh evbless fep_${k}; sleep 1; done

