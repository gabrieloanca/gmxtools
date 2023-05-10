#!/bin/bash

#ls fep_* | grep -v mdp | grep -v protocol | xargs rm 2> /dev/null
#ls restart* | grep -v mdp | xargs rm 2> /dev/null
#rm mdout.mdp run.log slurm* 2> /dev/null

#ls *.log | grep -v run.log | xargs rm 2> /dev/null
for j in rep_00?; do rm $j/*.edr; done
for j in rep_00?; do rm $j/fep*.cpt; done
for j in rep_00?; do rm $j/restart.cpt; done
for j in rep_00?; do rm $j/*.tpr; done
for j in rep_00?; do rm $j/fep*.log; done
for j in rep_00?; do rm $j/restart.log; done
for j in rep_00?; do rm $j/fep*.gro; done
for j in rep_00?; do rm $j/fep*.xvg; done
for j in rep_00?; do rm $j/restart.gro; done
for j in rep_00?; do rm $j/topol*; done
for j in rep_00?; do rm $j/mdout.mdp; done

for j in rep_00?; do rm $j/rerun/slurm*; done
for j in rep_00?; do rm $j/rerun/*.edr; done
for j in rep_00?; do rm $j/rerun/*.trr; done
for j in rep_00?; do rm $j/rerun/*.log; done
for j in rep_00?; do rm $j/rerun/fep*.xvg; done
for j in rep_00?; do rm $j/rerun/*.tpr; done
for j in rep_00?; do rm $j/rerun/mdout.mdp; done

