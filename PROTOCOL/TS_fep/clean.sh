#!/bin/bash

#ls fep_* | grep -v mdp | grep -v protocol | xargs rm 2> /dev/null
#ls restart* | grep -v mdp | xargs rm 2> /dev/null
#rm mdout.mdp run.log slurm* 2> /dev/null

ls *.log | grep -v run.log | xargs rm 2> /dev/null
rm *.edr
#rm fep*.cpt
rm *.tpr
rm fep*.log
#rm fep*.gro
rm topol*
rm mdout.mdp
rm fep*.xvg

rm rerun/slurm*
rm rerun/*.edr
rm rerun/*.trr
rm rerun/*.log
rm rerun/fep*.xvg
rm rerun/*.tpr
rm rerun/mdout.mdp

