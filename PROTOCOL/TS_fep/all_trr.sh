#!/bin/bash

export GMX_BIN=/home/x_gaboa/bin/gromacs-2022.4/bin/gmx

##gmx trjcat -f fep_000.trr fep_001.trr fep_002.trr fep_003.trr fep_004.trr fep_005.trr fep_006.trr fep_007.trr fep_008.trr fep_009.trr fep_010.trr fep_011.trr fep_012.trr fep_013.trr fep_014.trr fep_015.trr fep_016.trr fep_017.trr fep_018.trr fep_019.trr fep_020.trr fep_021.trr fep_022.trr fep_023.trr fep_024.trr fep_025.trr fep_026.trr fep_027.trr fep_028.trr fep_029.trr fep_030.trr fep_031.trr fep_032.trr fep_033.trr fep_034.trr fep_035.trr fep_036.trr fep_037.trr fep_038.trr fep_039.trr fep_040.trr fep_041.trr fep_042.trr fep_043.trr fep_044.trr fep_045.trr fep_046.trr fep_047.trr fep_048.trr fep_049.trr fep_050.trr -o all.trr -cat

## manual-2021, pg. 163
## -sort  : sorts trajectory files, not frames
## gmx trjcat concatenates several input trajectory files in sorted order
## Using -cat, you can simply paste several files together without removal
## of frames with identical time stamps


#${GMX_BIN} trjcat -f fep_0??.trr -o all.trr -sort -cat
echo {20,20} | ${GMX_BIN} trjconv -s equil_11_300k_npt.gro -f all.trr -o all.xtc -pbc nojump -skip 10 -n index.ndx -center
#echo 0 | ${GMX_BIN} trjconv -s fep_000.gro -f all.trr -o all_wat.xtc -pbc atom -ur compact -skip 10 -n index.ndx
#rm all.trr


