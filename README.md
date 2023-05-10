Python scripts for running EVB simulations Gromacs.

For download and updates, vizit or clone:
    https://github.com/gabrieloanca/gmxtools.git
    git@github.com:gabrieloanca/gmxtools.git
For suggestions, reporting bugs or for any assistance write to oanca.gabriel@gmail.com


The following scripts are available:
   ffld2gmx.py     - converts ffld parameters to OPLS-AA types fro Gromacs
   genposre.py     - generates posre files with different constraints for region 1 and region 2
   gmx2qfep.py     - writes outputs from 'gmx energy' into qfep5_gmx format
   gmx4evb.py      - builds topologies, one for each FEP frame
   qfep5_gmx       - Q5 mapping to extract the EVB profile (fortran code)
   qstats.py       - perfores statistics on qfep5_gmx output
   poly_data.py    - smoothen qfep5_gmx generated profiles by a 6th degree polynomial function

