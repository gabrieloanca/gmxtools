#!/bin/bash -l

#SBATCH --job-name=equil
#SBATCH -A snic2022-3-2
#SBATCH --time=2:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx
export RUNDIR=$TMPDIR/$SLURM_JOB_ID

echo "Host: $(hostname)"
echo "Tmpdir: $RUNDIR"
echo "Jobdir: $SLURM_SUBMIT_DIR"

prev='equil_11_300k_npt'
step='equil'
topol='topol_000.top'
#pdb='hbdh-ions.pdb'

# copy files to scratch dir
rsync ${step}.mdp ${RUNDIR}/
rsync ${prev}.cpt ${RUNDIR}/
rsync ${prev}.gro ${RUNDIR}/
rsync ${topol} ${RUNDIR}/
rsync ${pdb} ${RUNDIR}/ 
rsync -ah table* ${RUNDIR}/
rsync -r oplsaa.ff ${RUNDIR}/
rsync posre003.itp ${RUNDIR}/
cd ${RUNDIR}/

${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p topol_000.top -o ${step}.tpr -t ${prev}.cpt #-r ${pdb} #-maxwarn 1
time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16
rsync -ah $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"

## cleanup
#rm -rf $RUNDIR
echo "Done"

