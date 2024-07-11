#!/bin/bash -l

#SBATCH -A snic2022-3-2
#SBATCH --time=22:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx
export RUNDIR=$TMPDIR/$SLURM_JOB_ID

echo "Host: $(hostname)"
echo "Tmpdir: $RUNDIR"
echo "Jobdir: $SLURM_SUBMIT_DIR"

prev='fep_025'
step='ts_equil'
topol='topol_025.top'
pdb='tim-ions.pdb'

# copy files to scratch dir
#rsync ${step}.mdp ${RUNDIR}/
#rsync ${prev}.cpt ${RUNDIR}/
#rsync ${prev}.gro ${RUNDIR}/
#rsync ${topol} ${RUNDIR}/
#rsync ${pdb} ${RUNDIR}/ 
#rsync -ah table* ${RUNDIR}/
#rsync -r oplsaa.ff ${RUNDIR}/
#rsync posre_fep.itp ${RUNDIR}/
rsync -ah $SLURM_SUBMIT_DIR/ $RUNDIR/ --exclude="slurm*" --exclude="*.sh"
cd ${RUNDIR}/

${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${topol} -o ${step}.tpr -t ${prev}.cpt -r ${pdb} -maxwarn 1
time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16
rsync -ah $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"

## cleanup
#rm -rf $RUNDIR
echo "Done"

