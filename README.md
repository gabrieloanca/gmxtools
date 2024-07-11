**gmxtools** are a series of python scripts for EVB simulations with GROMACS.  
For download and updates, vizit or clone:  
>https://github.com/gabrieloanca/gmxtools.git or  
>git@github.com:gabrieloanca/gmxtools.git  

If you find these tools useful, please cite the following paper:  
> Gabriel Oanca, Florian van der Ent, Johan Åqvist, Efficient Empirical Valence Bond Simulations with GROMACS, *Journal of Chemical Theory and Computation*, **2023**, http://doi.org/10.1021/acs.jctc.3c00714  
  
For suggestions, reporting bugs or for any assistance write to gabriel.oanca@icm.uu.se or to oanca.gabriel@gmail.com.  

For easier use, add this folder to the **PATH** environment variable in the *bash rc* file.  
  
  
The following scripts are available:
>**ffld2gmx.py**      - converts ffld parameters to OPLS-AA types for Gromacs  
>**genposre.py**      - generates posre files with different constraints for the EVB and the non-EVB atoms  
>**gmx4evb.py**       - builds topologies, one for each FEP frame  
>**mapevb.py**        - analyses the energies and returns the EVB profile  
>**stats.py**         - calculates the mean and standard deviation over several replicas, from mapevb.py output files 
>**poly.py**          - smoothens the EVB profiles by a 6th degree polynomial  

To show this list in your terminal, type `gmxtools`.  
To get help for any of these tools, type `-h` after a tool's name (e.g., `gmx4evb.py -h`).  


• For how to write the *qmatoms.dat* file, follow the instruction inside *examples/HOW_TO_QMATOMS* file. The file *qmatoms.dat* can also be found inside the *examples/* folder.  
  
• A step-by-step guide for an EVB simulation in GROMACS can be found inside *examples/HOW_TO* file as well as in the Supporting Information of the paper above.  

### UPDATES APR. 11, 2024  
NOTE: These updates will slightly change the protocol and the workflow presented in the Supporting Information of the paper mentioned above. If you want to follow the workflow precisely as there, then download the tools and the protocol from **jctc_2023** branch (or check the **JCTC2023** tag). The old protocol can also be found inside **PROTOCOL(JCTC2023)** folder. The new changes are as follows:

 1. A new tool for mapping is introduced, **mapevb.py**, that will replace the **QFEP** tool of **Q** software (step III.4 in SI). By using **mapevb.py**, you can skip using **gmx2qfep.py** as well (step III.3 in SI).  
 2. **get_ene.sh** writes all the energy files in one folder named **ENE** (step. III.2). Execute it like `get_ene.sh 15 13`, where in this case '15' and '13' correspond to the numbers that will return the potential when calling the **energy** tool of GROMACS onto the files rerun with the first/last topolgy and with tge evbless.top, respectively.  
 3. The protocol changed so that the rerun is now performed after each FEP frame to avoid the requirement for large memory storage. The new protocol can be found inside the **PROTOCOL** folder (the old one is also saved inside **PROTOCOL(JCTC2023)** folder).  
 2. **get_ene.sh** writes all the energy files in one folder named **ENE** (step. III.2). Execute it like `get_ene.sh 15 13`, where in this case '15' and '13' correspond to the numbers that will return the potential when calling the **energy** tool of GROMACS onto the files rerun with the first/last topolgy and with evbless.top, respectively.  
 3. The protocol changed so that the rerun is now performed after each FEP frame to avoid the requirement for large memory storage by deleting the .trr files after the rerun. The new protocol can be found inside the **PROTOCOL** folder (the old one is also saved inside **PROTOCOL(JCTC2023)** folder).  
 4. **gmx4evb.py** now generates also the tabulated potential files together with the topologies, so you can skip the step I.17 in SI. When generating the tables, you can choose between single or double precision depending on your GROMACS installation (default: single), and the cut-off (default: 10Å).  
 5. For the same pair of atoms, you can now use soft-repulsions with different *$\beta$* values, one for each EVB state.  
    NOTE: For donor-acceptor pairs, you cannot combine soft-repulsion in one state with Lennard Jones in the other state.  
 6. The *.ac* file is now optional in **ffld2gmx.py**. If `--resp no` option is given or if the *.ac* file is missing, then *ffld_server* charges are used (step I.2 in SI). This option is not so important since charges must be written inside the *qmatoms.dat* file anyway, but it was annoying not to be able to skip the Gaussian calculations.  
