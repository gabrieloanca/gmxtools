#!/bin/bash

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

str='ts_equil.gro'

${GMX_BIN} grompp -f sysA.mdp    -c ${str} -p topol_000.top -o sysA.tpr    -maxwarn 2
${GMX_BIN} grompp -f sysB.mdp    -c ${str} -p topol_050.top -o sysB.tpr    -maxwarn 2
${GMX_BIN} grompp -f evbless.mdp -c ${str} -p evbless.top   -o evbless.tpr -maxwarn 2


