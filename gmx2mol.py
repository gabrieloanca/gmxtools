#!/usr/bin/env python3
# coding: utf-8

# #### Last update: June 27 2023
# #### It writes Gromacs generated data for an EVB job in gap
# #### files for the mapping tool of Molaris (mapping_hpc9.15).

# Beer-ware licence
# oanca.gabriel@gmail.com

# For download and updates, vizit or clone:
#     https://github.com/gabrieloanca/gmxtools.git  
#     git@github.com:gabrieloanca/gmxtools.git
# For suggestions, reporting buggs or for any assistance write to oanca.gabriel@gmail.com


import os, sys, re
import numpy as np
import getopt

def get_args():
    arg_list = sys.argv[1:]
    options="hf:r:p:o:s:"
    long_options=["help","frames=","reactants=","products=","outdir=","skip="]

    opts, args= getopt.getopt(arg_list, options, long_options)
    rs, ps, outdir, skip = 'sysA', 'sysB', 'gmx2mol', 1
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            show_help()
            sys.exit()
        elif opt in ('-f', '--frames'):
            frames = int(arg)
        elif opt in ('-r', '--reactants'):
            rs = str(arg)
        elif opt in ('-p', '--products'):
            ps = str(arg)
        elif opt in ('-o', '--outdir'):
            outdir = str(arg)
        elif opt in ('-s', '--skip'):
            skip = int(arg)
                
    return frames, rs, ps, outdir, skip


### for jupyter-notebook
### Set variables
##in_dir = 'fep_01/'   # for mol2mol
##out_dir = 'gmx2mol'
##rs = 'runA'
##ps = 'runB'
##frames = 51

##### READS MOLARIS GAP FILES
##def mol_gap_list(in_dir):
##    gaps = []
##    files = os.listdir(in_dir)   # remove "/fep_01" in map.py script
##    
##    for file in files:
##        if "gap" in file:
##            gaps.append(file)
##    
##    return sorted(gaps)###def read_mol_gaps(in_dir):
###    gaps = mol_gap_list(in_dir)
###    
###    # read the gap files and store the data
###    for i, gap in enumerate(gaps):
###        # set holders
###        by_line = []
###        h11_lines = []
###        h22_lines = []
###    
###        with open(f"{in_dir}/{gap}") as f:
###            data = f.read().strip().split("\n")
###    
###        # remove first entry and get lambda(RS)
###        data.pop(0)
###        
###        for line in data:
###            by_line.append(line.split())
###
###        for no in range(len(by_line)):
###            if not no%5:
###                h11_lines.append(by_line[no])
###                h22_lines.append(by_line[no+2])
###    
###        if i == 0:
###            npts = len(h11_lines)
###            frames = len(gaps)
###            rs_data = np.zeros(shape=(frames,npts))
###            ps_data = np.zeros(shape=(frames,npts))
###
###        # extract energy points
###        for j in range(npts):
###            rs_data[i][j] = float(h11_lines[j][0])
###            ps_data[i][j] = float(h22_lines[j][1])
###
###    return rs_data, ps_data
##### ENDS READINS MOLARIS GAP FILES


def get_feps(state, frames):
    gaps = []
    files = os.listdir(state)
    
    for file in files:
        # .+ means more that one unknown character
        # files are names as runA_000_sysA.xvg of runA_000_evbless.xvg
        if re.search("fep_.+_sys._pot.xvg", file) or re.search("fep_.+_evbless_pot.xvg", file):
            gaps.append(file)

    gaps = sorted(gaps)

    feps = {}
    for i in range (frames):
        for gap in gaps:
            no = '{:0>3}'.format(i)
            if f'fep_{no}_' in gap:
                if f'fep_{no}' in feps.keys():
                    feps[f"fep_{no}"].append(gap)
                else:
                    feps[f"fep_{no}"] = [gap]

    return feps


def read_file(file):
    lowlim = 24 # 24 first line are info
    with open(file) as f:
        data = f.read().strip().split("\n")[lowlim:]
        return data


# reads Gromacs generated data
def get_ene(feps, frames, state):
    heads = sorted(feps.keys())
    file = os.path.join(state, feps[heads[0]][0])
    rows = len(read_file(file))
    data = np.zeros(shape=(frames,rows))
    
    for i, head in enumerate(heads[:frames]):
        temp = np.zeros(shape=(len(feps[head]),rows))
        
        for j, file in enumerate(feps[head]):
            file = os.path.join(state, file)
            lst = read_file(file)
            if 'evbless' in file:
                for k,line in enumerate(lst):
                    temp[j][k] = - float(line.split()[1])
            else:
                for k,line in enumerate(lst):
                    temp[j][k] = float(line.split()[1])
        
        data[i]=np.sum(temp,axis=0)
        del(temp)
    return data * 0.239006


def write2gap(rs_data, ps_data, outdir, skip):
    gaps, npts = rs_data.shape
    try:
        os.mkdir(outdir) #create output directory
    except:
        pass
    dl = 1/(gaps-1)
    odir = os.path.abspath(outdir)
    
    for i in range(gaps):
        with open(f"{outdir}/map_evb.gap{i+1:0>3}", "w") as f:
            for j in range(npts):
                if not (j%skip):
                    rs = rs_data[i][j]
                    ps = ps_data[i][j]
                    l1 = 1 - i * dl     # lambda 1
                    l2 = i * dl         # lambda 2
                    pt = j * 10         # the point are recorded each 10 steps
                    if (i == gaps-1) and (j == 0):
                        f.write(f"""          {pt:.0f}   0.001000   0.00   2   {l1:.3f}   {l2:.3f}   0   {l2:.3f}   {dl:.4f}   {odir}/evb_lra.out   
   0.00   0.00   {rs:.2f}   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00
   0.00   0.00   0.00 0  0
   0.00   0.00   {ps:.2f}   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00
   0.00   0.00   0.00 0  0\n""")
                    elif j == npts - 1:
                        f.write(f"""          {pt:.0f}   0.001000   0.00   2   {l1:.3f}   {l2:.3f}   0   {l2:.3f}   {dl:.4f}
   0.00   0.00   {rs:.2f}   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00
   0.00   0.00   0.00 0  0
   0.00   0.00   {ps:.2f}   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00
   0.00   0.00   0.00 0  0""")                    
                    else:
                        f.write(f"""          {pt:.0f}   0.001000   0.00   2   {l1:.3f}   {l2:.3f}   0   {l2:.3f}   {dl:.4f}
   0.00   0.00   {rs:.2f}   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00
   0.00   0.00   0.00 0  0
   0.00   0.00   {ps:.2f}   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00\
   0.00   0.00   0.00   0.00
   0.00   0.00   0.00 0  0\n""")

###rs_data, ps_data = read_mol_gaps(indir)  # Molaris

def show_help():
    print("""
usage: gmx2mol.py [-h] -f #frames [-r RS folder] [-p PS folder] [-o output folder] [--skip skip]

gmx2mol.py takes the output of an EVB simulation in Gromacs and writes the energies for each FEP frame
into 'gap' files suited for being analyzed by mapping_hpc9.15 tool of Molaris software.

For info about how to analyze data in Molaris, see Molaris manual.

Options:
 -h    --help                show this help and exit
 -f    --frames    <int>     numbers of FEP frames to be analyzed
 -r    --reactants [sysA]    the folder name containing the files of RS state
 -p    --products  [sysB]    the folder name containing the files of PS state
 -o    --outdir    [gmx2mol] output folder containing the .gap files
 -s    --skip      [1]       collect data every 'skip' points
""")

if __name__ == "__main__":
    try:
        frames, rs, ps, outdir, skip = get_args()
    except:
        if (sys.argv[1:]) and (sys.argv[1:][0] in ('-h', '--help')):
            pass
        else:
            print('gmx2mol.py -f #frames [-r RS folder] [-p PS folder] [-o output folder] [--skip skip]')
            print('For help type "gmx2mol.py -h"')
        sys.exit()
    
    for i, state in enumerate((rs, ps)):
        feps = get_feps(state, frames)
        
        if i == 0:
            rs_data = get_ene(feps, frames, state)
        elif i == 1:
            ps_data = get_ene(feps, frames, state)
        
    write2gap(rs_data, ps_data, outdir, skip)

