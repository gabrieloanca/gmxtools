#!/bin/bash -l

#SBATCH -A snic2022-3-2
#SBATCH --time=02:00:00
#SBATCH -n 1
#SBATCH -c 8

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

export RUNDIR=$TMPDIR/$SLURM_JOB_ID
echo "Host: $(hostname)"
echo "Tmpdir: $RUNDIR"
echo "Jobdir: $SLURM_SUBMIT_DIR"

# copy files to scratch dir
rsync -ah $SLURM_SUBMIT_DIR/ $RUNDIR/ --exclude="slurm*" --exclude="*.sh"

mdfile='evbless'                     # sysA, sysB, evbless
topol='evbless.top'                  # topol_000.top, topol_050.top, evbless.top
structure='equil_fep31.gro'

cd $RUNDIR

echo "Time: $(date)" > ${mdfile}.txt
start_time=$(date +%s)
echo "I'm working in this directory:" >> ${mdfile}.txt
echo "   $PWD" >> ${mdfile}.txt

echo >> ${mdfile}.txt
echo "Steps:" >> ${mdfile}.txt
steps=$(ls -1v fep_???.trr | sed 's/.trr//')  
echo "${steps[@]}" >> ${mdfile}.txt
echo >> ${mdfile}.txt

if [[ ${mdfile} != 'evbless' ]]; then
  ${GMX_BIN} grompp -f ${mdfile}.mdp -c ${structure} -p ${topol} -o ${mdfile}.tpr #-maxwarn 1

  for step in ${steps[@]}; do
    echo -n "Time: $(date)     " >> ${mdfile}.txt
    echo -n "Running the step: ${step}    --    " >> ${mdfile}.txt
    if ${GMX_BIN} mdrun -s ${mdfile}.tpr -rerun ${step}.trr -table table.xvg -tableb table_b0.xvg table_b1.xvg -e ${step}_${mdfile}.edr -g ${step}_${mdfile}.log -o ${step}_${mdfile}.trr -dhdl ${step}_${mdfile}.xvg -pin on -ntmpi 1 -ntomp 8
      then echo "OK" >> ${mdfile}.txt
      rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
      else
      echo "FAILED" >> ${mdfile}.txt
      echo >> ${mdfile}.txt
      echo "Time: $(date)" >> ${mdfile}.txt
      exit 1
    fi
  done

elif [[ ${mdfile} == 'evbless' ]]; then
  ${GMX_BIN} grompp -f ${mdfile}.mdp -c ${structure} -p ${topol} -o ${mdfile}.tpr -maxwarn 1

  for step in ${steps[@]}; do
    echo -n "Time: $(date)     " >> ${mdfile}.txt
    echo -n "Running the step: ${step}    --    " >> ${mdfile}.txt
    if ${GMX_BIN} mdrun -s ${mdfile}.tpr -rerun ${step}.trr -e ${step}_${mdfile}.edr -g ${step}_${mdfile}.log -o ${step}_${mdfile}.trr -dhdl ${step}_${mdfile}.xvg -pin on -ntmpi 1 -ntomp 8
      then echo "OK" >> ${mdfile}.txt
      rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
    else
      echo "FAILED" >> ${mdfile}.txt
      echo >> ${mdfile}.txt
      echo "Time: $(date)" >> ${mdfile}.txt
      exit 1
    fi
  done
fi

echo >> ${mdfile}.txt
echo "Normal termination" >> ${mdfile}.txt
echo "Time: $(date)" >> ${mdfile}.txt
stop_time=$(date +%s)
echo "Total time: $((stop_time-start_time)) seconds" >> ${mdfile}.txt
#tt=$(bc <<< ${stop_time}-${start_time})
#h=$(bc <<< ${tt}/3600)
#rest=$(bc <<< ${tt}%3600)
#m=$(bc <<< ${rest}/60)
#s=$(bc <<< ${rest}%60)
#echo "Total time: ${h}:${m}:${s}" >> ${mdfile}.txt

rsync -ah --update $RUNDIR/ $SLURM_SUBMIT_DIR/ --exclude="*tmp"
## cleanup
rm -rf $RUNDIR
echo "Done"


