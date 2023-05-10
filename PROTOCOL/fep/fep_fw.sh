#!/bin/bash -l

#SBATCH --job-name=fep_fw
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

echo "Time: $(date)" > run_fw.log
start_time=$(date +%s)
echo -n "I'm working in this directory:" >> run_fw.log
echo "   $PWD" >> run_fw.log

echo >> run_fw.log
echo "Steps:" >>  run_fw.log				     
steps=$(ls -1v fep_{025..050}.mdp | sed 's/.mdp//')  
echo "${steps}" >> run_fw.log
echo >> run_fw.log
start_str='restart'
restraints='tim-ions.pdb'
prev='restart'

for j in {025..50}; do
  step=fep_${j}
  top=topol_${j}
  echo -n "Time: $(date)     " >> run_fw.log
  echo -n "Running the step: ${step}    --    " >> run_fw.log
  ${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${top}.top -o ${step}.tpr -t ${prev}.cpt -r ${restraints} -maxwarn 1
  if time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16; then
    echo "OK" >> run_fw.log
    rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
  else
    echo "FAILED" >> run_fw.log
    echo >> run_fw.log
    echo "Time: $(date)" >> run_fw.log
    break
  fi
  prev=${step}
done

echo >> run_fw.log
echo "Normal termination" >> run_fw.log
echo "Time: $(date)" >> run_fw.log
stop_time=$(date +%s)
echo "Total time: $((stop_time-start_time)) seconds" >> run_fw.log

rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
## cleanup
rm -rf $RUNDIR
echo "Done"

