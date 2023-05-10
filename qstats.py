#!/usr/bin/env python3
# coding: utf-8

# #### Date: Apr 11 2023
# #### This script extrats the EVB profiles from qfep output files (Part 3) and calculates the mean and std for dG# and dG0

# beer-ware licence
# oanca.gabriel@gmail.com

# For download and updates, vizit or clone:
#     https://github.com/gabrieloanca/gmxtools.git
#     git@github.com:gabrieloanca/gmxtools.git
# For suggestions, reporting bugs or for any assistance write to oanca.gabriel@gmail.com


import sys, os, re
import statistics as st
import matplotlib.pyplot as plt
import numpy as np

def names(fname):
    base = []
    pre = []
    
    try:
        pre = fname[:fname.index('.')]
    except:
        pre = fname
    
    for j, i in enumerate(pre):
        try:
            if type (int(i)) == int:
                pre[j] = '.'
        except:
            pass

    base = ''.join(pre)
    return base
    
def get_files(base):
    files = os.listdir()
    fdata=[]
    for f in files:
        if re.search(base, f):
            if (f[0] == '#'):
                continue
            elif 'graph' in f[f.index('.'):]:
                continue
            else:
                fdata.append(f)
    return fdata

def get_profiles(fdata):
    profiles = {}
    for i,profile in enumerate(fdata):
        profiles[profile] = {'x':[], 'y':[]}
        
        with open(profile) as f:
            data = f.read().strip().split("\n")
        
        for i, line in enumerate(data):
            if "# Part 3:" in line:
                click = i

        for line in data[click+2:]:
            l = line.strip().split()
            profiles[profile]['x'].append(float(l[1]))
            profiles[profile]['y'].append(float(l[3]))
    
        del(data)
    return profiles

def stats(profiles):
    ts_data, ps_data = [], []
    keys = list(profiles.keys())
    ts_loc = int(len(profiles[keys[0]]['x'])/2)
    for key in keys:
        rs = min(profiles[key]['y'][:ts_loc])
        rs_loc = profiles[key]['y'][:ts_loc].index(rs)
        ps = min(profiles[key]['y'][ts_loc:])
        ps_loc = profiles[key]['y'][ts_loc:].index(ps) + ts_loc
        ts = max(profiles[key]['y'][rs_loc:ps_loc])

        ps_data.append(ps)
        ts_data.append(ts)

    ts_mean = st.mean(ts_data)
    ts_std = st.stdev(ts_data)
    ps_mean = st.mean(ps_data)
    ps_std = st.stdev(ps_data)

    return ts_mean, ts_std, ps_mean, ps_std

def savefig(profiles, name='evb_profile'):
    keys = list(profiles.keys())
    fig = plt.figure()
    fig.set_figheight(3)
    fig.set_figwidth(4)
    ax1 = fig.add_axes([0, 0, 1, 1])
    ax1.grid()
    ax1.set_xlabel(r'$\epsilon_1$ - $\epsilon_2$ (kcal/mol)')
    ax1.set_ylabel(r'$\Delta$G (kcal/mol)')
    
    for key in keys:
        ax1.plot(profiles[key]['x'], profiles[key]['y'])
    fig.savefig(name,dpi = 300,orientation='landscape',facecolor='w',edgecolor='w',bbox_inches='tight')

def savedata(profiles):
    keys = list(profiles.keys())
    keys.sort()
    ts_loc = int(len(profiles[keys[0]]['x'])/2)
    #file = open(f"{name}.dat", "w")

    for key in keys:
        try:
            pre = key[:key.index('.')]
        except:
            pre = key
        file = open(f"{pre}.graph", "w")
        for i in range(len(profiles[key]['y'])):
            file.write(f"{profiles[key]['x'][i]}\t{profiles[key]['y'][i]:.2f}\n")
        file.close()


if __name__ == "__main__":
    try:
        fname = list(sys.argv[1])
    except:
        print()
        print('Usage: get_profiles.sh <data_file>')
        print('where "data_file" is any of the EVB data files.')
        print("Requirement: the names of the data files be numbered and must start with an alphabetic character. E.g.: rep0.dat, rep1.dat, etc.")
        sys.exit()

    base = names(fname)

    try:
        files = get_files(base)
    except:
        print(f'No files starting with "{base}" could be found.')
        print('Make sure that you are in the right directory and try again.')
        sys.exit()

    try:
        profiles = get_profiles(files)
    except:
        print("Data could not be extracted.")
        print("Make sure they have the right format and try again.")
        sys.exit()

    try:
        tmean, tstd, pmean, pstd = stats(profiles)

        print(f'\ndG#: {tmean:5.2f}  +/-{tstd:5.2f}')
        print(f'dG0: {pmean:5.2f}  +/-{pstd:5.2f}')
    except:
        print("\nEVB profiles could not be analyzed.")

    name='evb_profile.png'
    savefig(profiles, name)
    savedata(profiles)
    print(f'For a graphical representation, see "{name}".\n')



