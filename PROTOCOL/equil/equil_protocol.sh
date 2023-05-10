#!/bin/bash -l

#SBATCH --job-name=equil
#SBATCH -A snic2022-3-2
#SBATCH --time=20:00:00
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
echo "I'm working in this directory:" >> run.log
echo "   $PWD" >> run.log

echo >> run.log
echo "Steps:" >> run.log				     
steps=$(ls -1v equil*.mdp | sed 's/.mdp//')  
echo "${steps[@]}" >> run.log
pdb='tim-ions'
prev=''
echo >> run.log

for step in ${steps[@]}; do
    echo -n "Time: $(date)     " >> run.log
    echo -n "Running the step: ${step}    --    " >> run.log
    if [[ ${step} == "equil_01_5k_nvt" ]]; then
      ${GMX_BIN} grompp -f ${step}.mdp -c ${pdb}.pdb -p topol_000.top -o ${step}.tpr -r ${pdb}.pdb #-maxwarn 1
      if time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16 ; then
        echo "OK" >> run.log
        rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
      else 
        echo "FAILED" >> run.log
        echo >> run.log
        echo "Time: $(date)" >> run.log
        break
      fi
    elif [[ ${step} != "equil_01_5k_nvt" ]]; then
      ${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p topol_000.top -o ${step}.tpr -r ${pdb}.pdb -t ${prev}.cpt #-maxwarn 1
      if time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16 ; then
        echo "OK" >> run.log
        rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
      else
        echo "FAILED" >> run.log
        echo >> run.log
        echo "Time: $(date)" >> run.log
        break
      fi
    fi
    prev=${step}
done

echo >> run.log
echo "Normal termination" >> run.log
echo "Time: $(date)" >> run.log
stop_time=$(date +%s)
echo "Total time: $((stop_time-start_time)) seconds" >> run.log

## cleanup
rm -rf $RUNDIR
echo "Done"

