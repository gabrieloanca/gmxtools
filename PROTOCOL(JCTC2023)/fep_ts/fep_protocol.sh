#!/bin/bash -l

#SBATCH --job-name=fep
#SBATCH -A snic2022-3-2
#SBATCH --time=02:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

export RUNDIR=$TMPDIR/$SLURM_JOB_ID
echo "Host: $(hostname)"
echo "Tmpdir: $RUNDIR"
echo "Jobdir: $SLURM_SUBMIT_DIR"

# copy files to scratch dir
rsync -ah $SLURM_SUBMIT_DIR/ $RUNDIR/ --exclude="slurm*" --exclude="*.sh"
cd $RUNDIR

echo "Time: $(date)" > run.log
start_time=$(date +%s)
echo -n "I'm working in this directory:" >> run.log
echo "   $PWD" >> run.log

echo >> run.log
echo "Steps:" >>  run.log				     
steps=$(ls -1v fep_*.mdp | sed 's/.mdp//')  
echo "${steps[@]}" >> run.log
echo >> run.log
equil='restart'
start_str='equil_11_300k_npt.gro'
restraints='tim-ions.pdb'
prev='restart'

echo -n "Time: $(date)     " >> run.log
echo -n "Running the step: restart    --    " >> run.log
${GMX_BIN} grompp -f ${equil}.mdp -c ${start_str} -p topol_000.top -o ${equil}.tpr -r ${restraints}
if time ${GMX_BIN} mdrun -deffnm ${equil} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${equil}.gro -pin on -ntmpi 1 -ntomp 16; then
  echo "OK" >> run.log
  rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
else
  echo "FAILED" >> run.log
  echo >> run.log
  echo "Time: $(date)" >> run.log
fi

for j in {000..050}; do
  step=fep_${j}
  top=topol_${j}
  echo -n "Time: $(date)     " >> run.log
  echo -n "Running the step: ${step}    --    " >> run.log
  ${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${top}.top -o ${step}.tpr -t ${prev}.cpt -maxwarn 1
  if time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16; then
    echo "OK" >> run.log
    rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
  else
    echo "FAILED" >> run.log
    echo >> run.log
    echo "Time: $(date)" >> run.log
    break
  fi
  prev=${step}
done

echo >> run.log
echo "Normal termination" >> run.log
echo "Time: $(date)" >> run.log
stop_time=$(date +%s)
echo "Total time: $((stop_time-start_time)) seconds" >> run.log

rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
## cleanup
rm -rf $RUNDIR
echo "Done"

