#!/bin/bash

#sysses='sysA sysB evbless'
nosys=$1
noless=$2

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

echo > getene.log
for sys in ${sysses[@]}; do
  steps=$(ls -1v fep_???_${sys}.edr | sed 's/.edr//')
  echo >> getene.log
  echo ${steps[@]} >> getene.log
  echo >> getene.log

  if [[ ${sys} != "evbless" ]]; then
    if [ -d ${sys} ]; then
      rm -r ${sys}
      echo "${sys} directory was deleted" >> getene.log
      echo >> getene.log
    fi
    mkdir ${sys}

    for step in ${steps[@]}; do
      echo "${step}, ${no}, Potential, ${step}_pot.xvg" >> getene.log
      echo ${nosys} | ${GMX_BIN} energy -f ${step}.edr -o ${sys}/${step}_pot.xvg
      wait $!
    done
  else
    for step in ${steps[@]}; do
      echo "${step}, ${no}, Potential, ${step}_pot.xvg" >> getene.log
      echo ${noless} | ${GMX_BIN} energy -f ${step}.edr -o sysA/${step}_pot.xvg
      wait $!
      echo ${noless} | ${GMX_BIN} energy -f ${step}.edr -o sysB/${step}_pot.xvg
      wait $!
    done
  fi
done

# wait $!

