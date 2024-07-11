#!/bin/bash

sysses='sysA sysB evbless'
nosys=15
noless=13

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

if [ -d ENE ]; then
  rm -r ENE
fi
mkdir ENE

for sys in ${sysses[@]}; do
  steps=$(ls -1v fep_???_${sys}.edr | sed 's/.edr//')

  if [[ ${sys} != "evbless" ]]; then
    for step in ${steps[@]}; do
      echo ${nosys} | ${GMX_BIN} energy -f ${step}.edr -o ENE/${step}.xvg
      wait $!
    done
  else
    for step in ${steps[@]}; do
      echo ${noless} | ${GMX_BIN} energy -f ${step}.edr -o ENE/${step}.xvg
      wait $!
    done
  fi
done

