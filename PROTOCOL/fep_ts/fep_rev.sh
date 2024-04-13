#!/bin/bash -l

#SBATCH --job-name=fep_rev
#SBATCH -A snic2022-3-2
#SBATCH --time=01:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx
export RUNDIR=$TMPDIR/$SLURM_JOB_ID

echo "Host: $(hostname)"
echo "Tmpdir: $RUNDIR"
echo "Jobdir: $SLURM_SUBMIT_DIR"

# copy files to scratch dir
rsync -ah $SLURM_SUBMIT_DIR/ $RUNDIR/ --exclude="slurm*" --exclude="*.sh" --exclude="*.trr"
cd $RUNDIR

echo "Time: $(date)" > run_rev.log
start_time=$(date +%s)
echo -n "I'm working in this directory:" >> run_rev.log
echo "   $PWD" >> run_rev.log

echo >> run_rev.log
echo "Steps:" >>  run_rev.log				     
steps=$(ls -1v fep_{024..000}.mdp | sed 's/.mdp//')  
echo "${steps}" >> run_rev.log
echo >> run_rev.log
start_str=' restart'
restraints='tim-ions.pdb'
prev='restart'

for j in {024..000}; do
  step=fep_${j}
  top=topol_${j}
  echo -n "Time: $(date)     " >> run_rev.log
  echo -n "Running the step: ${step}    --    " >> run_rev.log
  ${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${top}.top -o ${step}.tpr -t ${prev}.cpt -r ${restraints} -maxwarn 1
  if time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16; then
    echo "OK" >> run_rev.log
    rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
  else
    echo "FAILED" >> run_rev.log
    echo >> run_rev.log
    echo "Time: $(date)" >> run_rev.log
    break
  fi
  prev=${step}
done

echo >> run_rev.log
echo "Normal termination" >> run_rev.log
echo "Time: $(date)" >> run_rev.log
stop_time=$(date +%s)
echo "Total time: $((stop_time-start_time)) seconds" >> run_rev.log

rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
## cleanup
rm -rf $RUNDIR
echo "Done"

