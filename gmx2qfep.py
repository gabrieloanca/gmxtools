#!/usr/bin/env python3
# coding: utf-8

# #### Last update: Dec 26 2022
# #### It writes Gromacs generated data for an EVB 
# #### job in files to be read by qfep tool of Q5.

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
    rs, ps, outdir, skip = 'sysA', 'sysB', 'qfep', 1
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            show_help()
            return 0, 0, 0, 0, 0
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


def get_feps(state, frames):
    gaps = []
    files = os.listdir(state)
    hlp = False
    
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
        l1 = 1 - i * dl     # lambda 1
        l2 = i * dl         # lambda 2
        with open(f"{outdir}/fep_{i:0>3}.dat", "w") as f:
            f.write(f"{l1:.6f}    {l2:.6f}\n")
            for j in range(npts):
                if not (j%skip):
                    rs = rs_data[i][j]
                    ps = ps_data[i][j]
                    f.write(f"{rs:.8f}    {ps:.8f}\n")
                    

def show_help():
    print("""
usage: gmx2qfep.py [-h] -f <#frames> [-r RS folder] [-p PS folder] [-o output folder] [--skip skip]

gmx2qfep.py takes the output of an EVB simulation in Gromacs and writes the energies into
files suited for being analyzed by qfep5_gmx tool. The qfep5_gmx tool is a modified version
of qfep5 tool of Q(5) software. For info about how to analyze the data with qfep, you can
follow the instruction from Q(5) manual (see: http://qdyn.no-ip.org/documents/qman.pdf)

Options:
 -h    --help             show this help and exit
 -f    --frames    <int>  numbers of FEP frames to be analyzed
 -r    --reactants [sysA] the folder name containing the files of RS state
 -p    --products  [sysB] the folder name containing the files of PS state
 -o    --outdir    [qfep] output folder containing the .gap files
 -s    --skip      [1]    collect data every <skip> points
""")


if __name__ == "__main__":
    try:
        frames, rs, ps, outdir, skip = get_args()
    except:
        print('gmx2qfep.py -f <#frames> [-r RS folder] [-p PS folder] [-o output folder] [--skip skip]')
        print('For help type "gmx2qfep.py -h"')
        sys.exit()
    
    if not([frames, rs, ps, outdir, skip] == [0,0,0,0,0]):
        for i, state in enumerate((rs, ps)):
            feps = get_feps(state, frames)
            
            if i == 0:
                rs_data = get_ene(feps, frames, state)
            elif i == 1:
                ps_data = get_ene(feps, frames, state)
            
        write2gap(rs_data, ps_data, outdir, skip)

