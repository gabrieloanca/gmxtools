#!/bin/bash -l

#SBATCH -A snic2022-3-2
#SBATCH --time=22:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

echo "Host: $(hostname)"
echo "Jobdir: $SLURM_SUBMIT_DIR"

prev='fep_025'
step='ts_equil'
topol='topol_025.top'
pdb='tim-ions.pdb'

${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${topol} -o ${step}.tpr -t ${prev}.cpt -r ${pdb} -n index.ndx -maxwarn 1
time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16

echo "Done"

