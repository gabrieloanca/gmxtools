#!/bin/bash

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

${GMX_BIN} trjcat -f fep_0??.xtc -o all_raw.xtc -sort -cat
#echo {20,20} | ${GMX_BIN} trjconv -s equil_11_300k_npt.gro -f all_raw.xtc -o all.xtc -pbc nojump -skip 10 -center #-n index.ndx
echo {20,20} | ${GMX_BIN} trjconv -s sysA.tpr -f all_raw.xtc -o all.xtc -center -ur compact -pbc mol -skip 10 #-n index.ndx

## manual-2021, pg. 163
## -sort  : sorts trajectory files, not frames
## gmx trjcat concatenates several input trajectory files in sorted order
## Using -cat, you can simply paste several files together without removal of frames with identical time stamps

