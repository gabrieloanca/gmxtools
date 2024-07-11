#!/bin/bash -l

#SBATCH --job-name=fep
#SBATCH -A snic2022-3-2
#SBATCH --time=02:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

log='run.log'
echo "Host: $(hostname)"
echo "Jobdir: $SLURM_SUBMIT_DIR"

echo "Time: $(date)" > ${log}
start_time=$(date +%s)
echo -n "I'm working in this directory:" >> ${log}
echo "   $PWD" >> ${log}
echo >> ${log}

equil='restart'
start_str='equil_11_300k_npt.gro'
restraints='tim-ions.pdb'
prev='restart'

echo -n "Time: $(date)     " >> ${log}
echo -n "Running the step: restart    --    " >> ${log}
${GMX_BIN} grompp -f ${equil}.mdp -c ${start_str} -p topol_000.top -o ${equil}.tpr -r ${restraints} -n index.ndx
if time ${GMX_BIN} mdrun -deffnm ${equil} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${equil}.gro -pin on -ntmpi 1 -ntomp 16; then
  echo "OK" >> ${log}
else
  echo "FAILED" >> ${log}
  echo >> ${log}
  echo "Time: $(date)" >> ${log}
fi

## RERUN tpr ##
${GMX_BIN} grompp -f sysA.mdp -c ${start_str} -p topol_000.top -o sysA.tpr -maxwarn 1
${GMX_BIN} grompp -f sysB.mdp -c ${start_str} -p topol_050.top -o sysB.tpr -maxwarn 1
${GMX_BIN} grompp -f evbless.mdp -c ${start_str} -p evbless.top -o evbless.tpr -maxwarn 1

for j in {000..050}; do
  step=fep_${j}
  top=topol_${j}
  echo -n "Time: $(date)     " >> ${log}
  echo -n "Running the step: ${step}    --    " >> ${log}
  ${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p ${top}.top -o ${step}.tpr -t ${prev}.cpt -n index.ndx -maxwarn 1
  if time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg table_b1.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16; then
    echo "OK" >> ${log}

    ## RERUN ##
    ${GMX_BIN} mdrun -s sysA.tpr -rerun ${step}.trr -table table.xvg -tableb table_b0.xvg table_b1.xvg -e ${step}_sysA.edr -g ${step}_sysA.log -o ${step}_sysA.trr -dhdl ${step}_sysA.xvg
    ${GMX_BIN} mdrun -s sysB.tpr -rerun ${step}.trr -table table.xvg -tableb table_b0.xvg table_b1.xvg -e ${step}_sysB.edr -g ${step}_sysB.log -o ${step}_sysB.trr -dhdl ${step}_sysB.xvg
    ${GMX_BIN} mdrun -s evbless.tpr -rerun ${step}.trr -e ${step}_evbless.edr -g ${step}_evbless.log -o ${step}_evbless.trr -dhdl ${step}_evbless.xvg
    rm ${step}.trr
  else
    echo "FAILED" >> ${log}
    echo >> ${log}
    echo "Time: $(date)" >> ${log}
    break
  fi
  prev=${step}
done

echo >> ${log}
echo "Normal termination" >> ${log}
echo "Time: $(date)" >> ${log}
stop_time=$(date +%s)
echo "Total time: $((stop_time-start_time)) seconds" >> ${log}
echo "Done"

