#!/usr/bin/env python3
# coding: utf-8

## Sept 1 2023

## This script reads QFEP output files and return dG# and dG0 plus STD from lambda profiles (under 'Part 1')

## Usage: ./qstats_lambda.py <input file>

import sys, os, re
import statistics as st
import matplotlib.pyplot as plt

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
            elif ('graph' in f[f.index('.'):]) or ('lambda' in f[f.index('.'):]):
                continue
            else:
                fdata.append(f)

    fdata.sort()
    return fdata

def get_profiles(fdata):
    profiles = {}
    for i,profile in enumerate(fdata):
        profiles[profile] = {'x':[], 'y':[]}
        
        with open(profile) as f:
            data = f.read().strip().split("\n")
        
        for i,l in enumerate(data):
            if "Part 1:" in l:
                start1 = i+4
            elif "Part 2:" in l:
                start2 = i+2
                stop1 = i-7

        for line in data[start1:stop1+1]:
            l = line.strip().split()
            profiles[profile]['x'].append(float(l[0]))
            profiles[profile]['y'].append(float(l[5]))
    
        del(data)
    return profiles

def stats(profiles):
    ts_data, ps_data = [], []
    keys = list(profiles.keys())
    for key in keys:
        ps = profiles[key]['y'][-1]
        ts = max(profiles[key]['y'])

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
    ax1.set_xlabel(r'$\lambda$')
    ax1.set_ylabel(r'$\Delta$G (kcal/mol)')
    
    for key in keys:
        ax1.plot(profiles[key]['x'], profiles[key]['y'])
    ax1.invert_xaxis()
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
        file = open(f"{pre}.lambda", "w")
        for i in range(len(profiles[key]['y'])):
            file.write(f"{profiles[key]['x'][i]: >7.3f}\t{profiles[key]['y'][i]: >8.3f}\n")
        file.close()

if __name__ == "__main__":
    try:
        fname = list(sys.argv[1])
    except:
        print()
        print('Usage: qstats_lambda.py <data_file>')
        print('where "data_file" is any of the EVB data files.')
        print("Requirement: the names of the data files be numbered and must start with an alphabetic character. E.g.: rep0.dat, rep1.dat, etc.")
        print()
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
        print("\nlambda profiles could not be analyzed.")
        sys.exit()

    name='lambda_profiles.png'
    savefig(profiles, name)
    savedata(profiles)
    print(f'For a graphical representation, see "{name}".\n')
