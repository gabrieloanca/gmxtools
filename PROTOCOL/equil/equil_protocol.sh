#!/bin/bash -l

#SBATCH --job-name=equil
#SBATCH -A snic2022-3-2
#SBATCH --time=20:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --gpus-per-task=1

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

log='run.log'
echo "Host: $(hostname)"
echo "Jobdir: $SLURM_SUBMIT_DIR"

echo "Time: $(date)" > ${log}
start_time=$(date +%s)
echo "I'm working in this directory:" >> ${log}
echo "   $PWD" >> ${log}

echo >> ${log}
echo "Steps:" >> ${log}
steps=$(ls -1v equil*.mdp | sed 's/.mdp//')  
echo "${steps[@]}" >> ${log}
echo >> ${log}

pdb='tim'
prev='minimization'

${GMX_BIN} grompp -f minimization.mdp -c ${pdb}.pdb -p topol_000.top -o minimization.tpr -r ${pdb}.pdb #-maxwarn 1
echo -n "Time: $(date)     " >> ${log}
echo -n "Running the step: minimization    --    " >> ${log}
if time ${GMX_BIN} mdrun -deffnm minimization -table table.xvg -tableb table_b0.xvg -c minimization.gro -pin on -ntmpi 1 -ntomp 16 ; then
    echo "OK" >> ${log}
else 
    echo "FAILED" >> ${log}
    echo >> ${log}
    echo "Time: $(date)" >> ${log}
    break
fi

for step in ${steps[@]}; do
    echo -n "Time: $(date)     " >> ${log}
    echo -n "Running the step: ${step}    --    " >> ${log}
    ${GMX_BIN} grompp -f ${step}.mdp -c ${prev}.gro -p topol_000.top -o ${step}.tpr -r ${prev}.gro -t ${prev}.cpt -n index.ndx #-maxwarn 1
    if time ${GMX_BIN} mdrun -deffnm ${step} -table table.xvg -tableb table_b0.xvg -c ${step}.gro -pin on -ntmpi 1 -ntomp 16 ; then
        echo "OK" >> ${log}
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

