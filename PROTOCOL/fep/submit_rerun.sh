#!/bin/bash

for i in rep_00?; do cd ${i}/rerun/; for k in {000..050}; do sbatch rerun.sh sysA fep_${k}; sleep 0.5; done; cd ../../; done

for i in rep_00?; do cd ${i}/rerun/; for k in {000..050}; do sbatch rerun.sh sysB fep_${k}; sleep 0.5; done; cd ../../; done

for i in rep_00?; do cd ${i}/rerun/; for k in {000..050}; do sbatch rerun.sh evbless fep_${k}; sleep 0.5; done; cd ../../; done

