#!/bin/bash -l

#SBATCH --job-name=rerun
#SBATCH -A snic2022-3-2
#SBATCH --time=00:20:00
#SBATCH -n 1
#SBATCH -c 4

mdfile=$1
step=$2

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx
export RUNDIR=${TMPDIR}/${SLURM_JOB_ID}

echo "Host: $(hostname)"
echo "Tmpdir: ${RUNDIR}"
echo "Jobdir: ${SLURM_SUBMIT_DIR}"
echo "Now running: ${mdfile}, ${step}"

# copy files to scratch dir
rsync ${SLURM_SUBMIT_DIR}/${mdfile}.tpr ${RUNDIR}/
rsync -a ${SLURM_SUBMIT_DIR}/${step}.trr ${RUNDIR}/
rsync ${SLURM_SUBMIT_DIR}/table* ${RUNDIR}/
rsync -r ${SLURM_SUBMIT_DIR}/oplsaa.ff ${RUNDIR}/

cd ${RUNDIR}/

if [[ ${mdfile} != 'evbless' ]]; then
    ${GMX_BIN} mdrun -s ${mdfile}.tpr -rerun ${step}.trr -table table.xvg -tableb table_b0.xvg -e ${step}_${mdfile}.edr -g ${step}_${mdfile}.log -o ${step}_${mdfile}.trr -dhdl ${step}_${mdfile}.xvg -pin on -ntmpi 1 -ntomp 4
elif [[ ${mdfile} == 'evbless' ]]; then
    ${GMX_BIN} mdrun -s ${mdfile}.tpr -rerun ${step}.trr -e ${step}_${mdfile}.edr -g ${step}_${mdfile}.log -o ${step}_${mdfile}.trr -dhdl ${step}_${mdfile}.xvg -pin on -ntmpi 1 -ntomp 4
fi

rsync -ah --update ${RUNDIR}/ ${SLURM_SUBMIT_DIR}/

## cleanup
rm -rf ${RUNDIR}
echo "Done"

