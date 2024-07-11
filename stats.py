#!/usr/bin/env python3
# coding: utf-8

# #### Last update: Apr 8 2024
# #### This script extrats the EVB profiles from mapevb.py output files and calculates the mean and std for dG# and dG0

# beer-ware licence
# oanca.gabriel@gmail.com

# For download and updates, vizit or clone:
#     https://github.com/gabrieloanca/gmxtools.git
#     git@github.com:gabrieloanca/gmxtools.git
# For suggestions, reporting bugs or for any assistance write to oanca.gabriel@gmail.com


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
            elif 'graph' in f[f.index('.'):]:
                continue
            elif 'inp' in f[f.index('.'):]:
                continue
            elif 'png' in f[f.index('.'):]:
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

        for line in data[1:]:
            l = line.strip().split()
            profiles[profile]['x'].append(float(l[0]))
            profiles[profile]['y'].append(float(l[3]))
    
        del(data)
    return profiles

def stats(profiles):
    ts_data, ps_data = [], []
    rs_list = {}  # I save these to shift the profiles since <dGg norm> in QFEP gets shifted at zero at the lowest state, which is not always RS
    keys = list(profiles.keys())
    ts_loc = int(len(profiles[keys[0]]['x'])/2)
    for key in keys:
        prev = profiles[key]['y'][0]
        for j, i in enumerate(profiles[key]['y'][1:ts_loc]):
            if i < prev:
                prev = i
            elif i > prev:
                rs = prev
                rs_loc = j + 1
                break

        #rs = min(profiles[key]['y'][:ts_loc])
        #rs_loc = profiles[key]['y'][:ts_loc].index(rs)

        prev = profiles[key]['y'][-1]
        for j, i in enumerate(profiles[key]['y'][-2:rs_loc:-1]):
            if i < prev:
                prev = i
            elif i > prev:
                ps = prev
                ps_loc = len(profiles[key]['y']) - j - 2
                break

        #ps = min(profiles[key]['y'][ts_loc:])
        #ps_loc = profiles[key]['y'][ts_loc:].index(ps) + ts_loc
        ts = max(profiles[key]['y'][rs_loc:ps_loc])

        rs_list[key] = rs
        ps_data.append(ps-rs)
        ts_data.append(ts-rs)

    ts_mean = st.mean(ts_data)
    try:
        ts_std = st.stdev(ts_data)
    except:
        ts_std = 0.0

    ps_mean = st.mean(ps_data)
    try:
        ps_std = st.stdev(ps_data)
    except:
        ps_std = 0.0

    return ts_mean, ts_std, ps_mean, ps_std, rs_list

def savefig(profiles, rs_list, name='evb_profile'):
    keys = list(profiles.keys())
    fig = plt.figure()
    fig.set_figheight(3)
    fig.set_figwidth(4)
    ax1 = fig.add_axes([0, 0, 1, 1])
    ax1.grid()
    ax1.set_xlabel(r'$\epsilon_1$ - $\epsilon_2$ (kcal/mol)')
    ax1.set_ylabel(r'$\Delta$G (kcal/mol)')
    
    for key in keys:
        ax1.plot(profiles[key]['x'], [i - rs_list[key] for i in profiles[key]['y']])
    fig.savefig(name,dpi = 300,orientation='landscape',facecolor='w',edgecolor='w',bbox_inches='tight')

def savedata(profiles, rs_list):
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
            file.write(f"{profiles[key]['x'][i]}\t{profiles[key]['y'][i] - rs_list[key]:.2f}\n")
        file.close()

def help():
    print('Usage: stats.py <data_file>')
    print('where "data_file" is any of the mappevb.py .dat files.')
    print("Requirements: the names of the data files be numbered and must start with an alphabetic character. E.g.: rep0.dat, rep1.dat, etc.")
    print()

if __name__ == "__main__":
    try:
        fname = list(sys.argv[1])
        if ''.join(fname) == '-h' or ''.join(fname) == '--help':
            sys.exit()
    except:
        print()
        help()
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
        tmean, tstd, pmean, pstd, rs_list = stats(profiles)
        print(f'\ndG#: {tmean:5.2f}  +/-{tstd:5.2f}')
        print(f'dG0: {pmean:5.2f}  +/-{pstd:5.2f}')
    except:
        print("\nEVB profiles could not be analyzed.")
        help()
        sys.exit()

    name='evb_profile.png'
    savefig(profiles, rs_list, name)
    #savedata(profiles, rs_list)
    print(f'For a graphical representation, see "{name}".\n')



