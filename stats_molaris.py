#!/usr/bin/env python3
# coding: utf-8

# #### Last_update: Dec 25 2022
# #### This scripts extrats dG# and dG0 from 'dG_dE.graph' of Molaris' mapping file

# beer-ware licence
# oanca.gabriel@gmail.com

# For download and updates, vizit or clone:
#     https://github.com/gabrieloanca/gmxtools.git
#     git@github.com:gabrieloanca/gmxtools.git
# For suggestions, reporting bugs or for any assistance write to oanca.gabriel@gmail.com


import sys, os
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
    
    for i in pre:
        try:
            if type (int(i)) == int:
                break
        except:
            base.append(i)    
    base = ''.join(base)
    return base
    
def get_files(base):
    files = os.listdir()
    fdata=[]
    for f in files:
        if base in f:
            fdata.append(f)
    return fdata

def get_profiles(fdata):
    profiles = {}
    for i,profile in enumerate(fdata):
        profiles[profile] = {'x':[], 'y':[]}
        
        with open(profile) as f:
            data = f.read().strip().split("\n")
        
        for l in data:
            l = l.strip().split()
            profiles[profile]['x'].append(float(l[0]))
            profiles[profile]['y'].append(float(l[1]))
    
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
        
        ps_data.append(ps - rs)
        ts_data.append(ts - rs)

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
    
    ts_loc = int(len(profiles[keys[0]]['x'])/2)
    for key in keys:
        rs = min(profiles[key]['y'][:ts_loc])
        y = np.array(profiles[key]['y']) - rs
        ax1.plot(profiles[key]['x'], y)
    fig.savefig(name,dpi = 300,orientation='landscape',facecolor='w',edgecolor='w',bbox_inches='tight')


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
        print("Make sure they have the right format of Molaris' mapping output and try again.")
        sys.exit()

    try:
        tmean, tstd, pmean, pstd = stats(profiles)

        print(f'\ndG#: {tmean:5.2f}  +/-{tstd:5.2f}')
        print(f'dG0: {pmean:5.2f}  +/-{pstd:5.2f}')
    except:
        print("\nEVB profiles could not be analyzed.")

    name='evb_profile.png'
    savefig(profiles, name)
    print(f'For a graphical representation, see "{name}".\n')



