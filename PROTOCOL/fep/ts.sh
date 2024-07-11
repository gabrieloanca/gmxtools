#!/bin/bash -l

#SBATCH --job-name=ts
#SBATCH -A snic2022-3-2
#SBATCH --time=00:20:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

step='restart'
prev='ts_equil'
topol='topol_025.top'
restraints='tim-ions.pdb'

## RERUN tpr ##
${GMX_BIN} grompp -f sysA.mdp -c ${prev}.gro -p topol_000.top -o sysA.tpr -maxwarn 1
${GMX_BIN} grompp -f sysB.mdp -c ${prev}.gro -p topol_050.top -o sysB.tpr -maxwarn 1
${GMX_BIN} grompp -f evbless.mdp -c ${prev}.gro -p evbless.top -o evbless.tpr -maxwarn 1
## END RERUN tpr ##

${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${topol} -o ${step}.tpr -r ${restraints} -n index.ndx -maxwarn 1
${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16

echo "Done"

