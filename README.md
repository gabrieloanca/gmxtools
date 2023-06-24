gmxtools is a series of python scripts for running EVB simulations in gromacs.

For download and updates, vizit or clone:
    https://github.com/gabrieloanca/gmxtools.git
    git@github.com:gabrieloanca/gmxtools.git
For suggestions, reporting bugs or for any assistance write to oanca.gabriel@gmail.com

For easier use, ddd the path to this folder into your bash rc file.


The following scripts are available:
   ffld2gmx.py     - converts ffld parameters to OPLS-AA types fro Gromacs
   genposre.py     - generates posre files with different constraints for region 1 and region 2
   gmx2qfep.py     - writes outputs from 'gmx energy' into qfep5_gmx format
   gmx4evb.py      - builds topologies, one for each FEP frame
   qfep5_gmx       - Q5 mapping to extract the EVB profile (fortran code)
   qstats.py       - extracts the mean and the standard deviation from qfep5_gmx output files
   poly_data.py    - smoothens qfep5_gmx EVB generated data by a 6th degree polynomial function


To show this list in your terminal, type *gmxtools*

To get help for any of these tools, type -h after a tool's name (e.g., gmx4evb.py -h)

