#!/bin/bash -l

#for j in *.mdp; do
	#sed -i "s/verlet-buffer-tolerance  = 0.005/verlet-buffer-tolerance  = 2e-04/" $j
	#sed -i "s/vdw-type                 = PME/vdw-type                 = cut-off/" $j
	#sed -i "s/fourierspacing           = 0.125/fourierspacing           = 0.12/" $j
	#sed -i "s/pme-order                = 4/pme-order                = 5/" $j
	#sed -i "s/coulomb-modifier         = Potential-shift ; shift el pot such that 0 at cut-off/coulomb-modifier         = None/" $j

	#sed "s/nstvout                  = 10000     ; saving velocities/nstvout                  = 0         ; saving velocities/" < $j > $j.tmp
	#sed "s/nstfout                  = 10000     ; saving forces/nstfout                  = 0         ; saving forces/" < $j > $j.tmp
	#sed "s/nstlog                   = 10000     ; write en to log file/nstlog                   = 0         ; write en to log file/" < $j > $j.tmp
	#sed "s/nstenergy                = 10000     ; fq to write en to .edr file/nstenergy                = 0         ; fq to write en to .edr file/" < $j > $j.tmp
	#sed "s/comm-mode                = none/comm-mode                = Linear/" < $j > $j.tmp
#	sed "s/nsteps                   = 10000/nsteps                   = 20000/" < $j > $j.tmp
	#sed -i "s///" $j
	#sed -i "s///" $j
	#sed -i "s///" $j
#done

#sed -i "s/init-lambda-state        = 23/init-lambda-state        = 22/" restart.mdp

for j in *.tmp; do
    mv $j `basename $j .tmp`
done

### sed -i "s///" $j
