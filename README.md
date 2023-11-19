**gmxtools** are a series of python scripts for EVB simulations with GROMACS.  
For download and updates, vizit or clone:  
>https://github.com/gabrieloanca/gmxtools.git or  
>git@github.com:gabrieloanca/gmxtools.git  

If you find these tools useful, please cite the following paper:  
> Gabriel Oanca, Florian van der Ent, Johan Åqvist, Efficient Empirical Valence Bond Simulations with GROMACS, *Journal of Chemical Theory and Computation*, **2023**, doi: 10.1021/acs.jctc.3c00714  

For suggestions, reporting bugs or for any assistance write to gabriel.oanca@icm.uu.se or to oanca.gabriel@gmail.com.  

For easier use, add this folder to the **PATH** environment variable in the bash rc file.  


The following scripts are available:  
>**ffld2gmx.py**  - converts ffld parameters to OPLS-AA types for Gromacs  
>**genposre.py**  - generates posre files with different constraints for region 1 and region 2  
>**gmx2qfep.py**  - writes the output from 'gmx energy' into qfep5_gmx format  
>**gmx4evb.py**   - builds topologies, one for each FEP frame  
>**qfep5_gmx**    - Q5 mapping to extract the EVB profile (fortran code)  
>**qstats.py**    - calculates the mean and the standard deviation from qfep5_gmx output files  
>**poly.py**      - smoothens qfep5_gmx EVB generated data with a 6th degree polynomial function  


To show this list in your terminal, type `gmxtools`.  
To get help for any of these tools, type `-h` after a tool's name (e.g., `gmx4evb.py -h`).  
  
  
• For how to write the *qmatoms.dat* file, follow the instruction inside *examples/HOW_TO_QMATOMS* file.
The file *qmatoms.dat* can also be found inside the *examples/* folder.

• A step-by-step guide for an EVB simulation in GROMACS can be found inside *examples/HOW_TO* file as well as in the Supporting Information of the paper above.
