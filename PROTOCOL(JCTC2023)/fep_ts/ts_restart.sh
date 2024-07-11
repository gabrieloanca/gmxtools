#!/bin/bash -l

#SBATCH --job-name=ts
#SBATCH -A snic2022-3-2
#SBATCH --time=00:20:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx
export RUNDIR=$TMPDIR/$SLURM_JOB_ID

step='restart'
prev='ts_equil'
topol='topol_025.top'
restraints='tim-ions.pdb'

# copy files to scratch dir
#rsync -ah $SLURM_SUBMIT_DIR/ $RUNDIR/ --exclude="slurm*" --exclude="*.sh"
rsync ${step}.mdp ${RUNDIR}/
#rsync ${prev}.cpt ${RUNDIR}/
rsync ${prev}.gro ${RUNDIR}/
rsync ${topol} ${RUNDIR}/
rsync ${restraints} ${RUNDIR}/ 
rsync -ah table* ${RUNDIR}/
rsync -r oplsaa.ff ${RUNDIR}/
rsync posre_fep.itp ${RUNDIR}/

cd $RUNDIR

${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${topol} -o ${step}.tpr -r ${restraints} -maxwarn 1
${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16
rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"

## cleanup
rm -rf $RUNDIR
echo "Done"

