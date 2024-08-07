#!/usr/bin/env python3
# coding: utf-8

# #### Last update: Apr 8 2023
# #### It creates a 6 degree polynomial profiles
# #### fitted onto the output data of mapevb.py

# beer-ware licence
# oanca.gabriel@gmail.com

# For download and updates, vizit or clone:
#     https://github.com/gabrieloanca/gmxtools.git  
#     git@github.com:gabrieloanca/gmxtools.git
# For suggestions, reporting bugs or for any assistance write to oanca.gabriel@gmail.com

import os, re
import argparse
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

def get_args():
    parser = argparse.ArgumentParser(epilog='''\
By default, it reads any file of the form rep_*.dat in the working directory.
If the files have different names, like fep_*.dat, then pass any of the input files to poly.py

poly.py -i fep_000 -o logfile.log
''')
    parser.add_argument("-s", "--left", help="skip point to the left (default = 3)", required=False, type=int, default=3)
    parser.add_argument("-r", "--right", help="skip point to the right (default = 3)", required=False, type=int, default=3)
    parser.add_argument("-i", "--input", help="Any of the input files. Default: 'rep_000'.", required=False, type=str, default='rep_000.dat')
    parser.add_argument("-o", "--output", help="Output name for the log file. Default: 'logfile.log'.", required=False, type=str, default='logfile.log')
    args = parser.parse_args()

    return args.left, args.right, args.input, args.output

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

def get_date():
    return datetime.now().strftime("%b/%d/%Y %H:%M:%S")

def get_files(base):
    files = os.listdir()
    fdata=[]
    for f in files:
        if re.search(base, f) and re.search(".dat", f):
            if (f[0] == '#') or (f[0] == '.'):
                continue
            else:
                fdata.append(f)
    
    fdata.sort()
    return fdata

def get_data(file):
    """
    read data from the output of mapevb.py and store it inside two numpy arrays for scikit-learn
    """
    try:
        with open(f"{file}") as f:
            data = f.read().split("\n")
    except:
        print(f"File {file} not found")

    # variables to store (e2-e1) and dG
    de, dg = [], []

    for l in data[1:]:
        l = l.strip().split()
        if l and l[0][0] != "#":
            i,j = l[0], l[-1]
            de.append(float(i))
            dg.append(float(j))
        else:
            pass
        
    de=np.array(de)
    dg=np.array(dg)
    
    return de, dg

def fit_poly(de, dg):
    """
    takes in the mapevb.py data and generates a 6 degree polinomial curve
    whose coefficients are determined by linear regression
    """
    de = de.reshape(-1, 1) 
    dg = dg.reshape(-1, 1)

    poly = PolynomialFeatures(degree = 6)
    x_poly = poly.fit_transform(de)
    lin = LinearRegression()
    lin.fit(x_poly, dg)

    predict = lin.predict(x_poly)
    return predict, lin.coef_

## bilds more fitted points than in the original data
#def fine_poly(de, coef):
#    #new_x = np.arange(de[0],de[-1],1)
#    new_x = np.linspace(de[0], de[-1], 270)
#    new_y = list(map(lambda x: coef[0][0] + x*coef[0][1] + x**2*coef[0][2] +\
#                               x**3*coef[0][3] + x**4*coef[0][4] + x**5*coef[0][5] +\
#                               x**6*coef[0][6], new_x))
#  
#    new_x = np.array(new_x)
#    new_y = np.array(new_y) 
#    return new_x, new_y

def stats(dg0s, dgts):
    dg0s = np.array(dg0s)
    dgts = np.array(dgts)
    dg_zero = dg0s.mean()
    std_zero = dg0s.std()
    dg_dagger = dgts.mean()
    std_dagger = dgts.std()
    
    return dg_zero, std_zero, dg_dagger, std_dagger

def save_report(dg_dagger, std_dagger, dg_zero, std_zero, reps, out): 
    """
    Saves std and mean for activation and reaction free energies
    """
    date = get_date()
    with open(out, "w") as file:
        file.write(f"""\n Date: {date}

 The following replicas have been analyzed:\n""")
        for rep in reps:
            file.write('    ' + rep + '\n')
                   
        file.write(f"""
   ΔG‡: {dg_dagger:.2f} ± {std_dagger:.2f}    SEM = {std_dagger/np.sqrt(len(reps)):>4.2f}
   ΔG°: {dg_zero:.2f} ± {std_zero:.2f}    SEM = {std_zero/np.sqrt(len(reps)):>4.2f}
=====================================""")

class Fig():
    """
    Plot more profiles, added in separate steps, on the same canvas
    """
    def __init__(self, ylim = ()):
        self.ylim = ylim
        self.fig = plt.figure()
        self.ax1 = self.fig.add_axes([0, 0, 1, 1])
        self.ax1.set_xlabel('ε₁-ε₂ (kcal/mol)')
        self.ax1.set_ylabel('ΔG (kcal/mol)')
        self.ax1.grid()
        if self.ylim == ():
            pass
        else:
            (self.low, self.high) = self.ylim
            self.ax1.set_ylim(self.low, self.high)

    # plot_dg() plots only the mapevb.py data
    def plot(self, x=[], y=[]):
        self.ax1.plot(x, y)

    # plot_poly() plots the polynimial fitted profile
    def plot_show(self):
        self.fig.show()

    # print_fig() saves the plot on a .png file
    def print_fig(self,name = 'profile.png'):
        #self.fig.savefig(f"{name}.png", dpi = 300, orientation = 'landscape',
        #                 format = 'png', transparent = True) 
        self.fig.savefig(f"{name}", dpi = 300, orientation = 'landscape',
                         facecolor = 'w', edgecolor = 'w', bbox_inches = 'tight')

if __name__ == "__main__":
    sx, dx, inp, out = get_args()

    base = names(inp)

    files = get_files(base)
    fig = Fig()

    dg0s, dgts = [], [] # dg0 and dg#
    chk = True
    
    for file in files:
        de, dg = get_data(file)  
        fit, coef = fit_poly(de, dg)
        fit = fit.reshape(1, -1)[0]
        #x, y = fine_poly(de, coef)
        mins, maxs = [], [] # minima and maxima on the graph
        
        if dx > 0:
            x = de[sx:-dx]
            y = fit[sx:-dx]
        else:
            x = de[sx:]
            y = fit[sx:]
            
        dydx = np.diff(y)/np.diff(x)
        k = dydx[0]
        for i,j in enumerate(dydx[1:]):
            if k<=0 and j > 0:
                mins.append((i+1,k))   # i+1 because i start from 1
            elif k>=0 and j< 0:
                maxs.append((i+1,k))
            k = j
        
        rs = 0
        if (len(mins) != 2) or (len(maxs) != 1):
            chk = False
            print(f"There are {len(mins)} minima and {len(maxs)} maxima in {file}.")
            print("Chekc fig.png graph and trim the data appropriately.")
            print()
        else:
            rs = y[mins[0][0]]
            dg0s.append(y[mins[1][0]]-rs)
            dgts.append((y[maxs[0][0]]-rs))

        fig.plot(x,y-rs)   
        #fig.plot_dg(x,y)

    if chk:
        dg_zero, std_zero, dg_dagger, std_dagger = stats(dg0s, dgts)
        print()
        print(f"   ΔG‡: {dg_dagger:>5.2f} ± {std_dagger:>4.2f}    SEM = {std_dagger/np.sqrt(len(files)):>3.2f}")
        print(f"   ΔG°: {dg_zero:>5.2f} ± {std_zero:>4.2f}    SEM = {std_zero/np.sqrt(len(files)):>3.2f}")

        if out == 'logfile.log':
            pass
        else:
            if '.' in out:
                pass
            else:
                out = out+'.log'

        save_report(dg_dagger, std_dagger, dg_zero, std_zero, files, out)

    if 'logfile' in out: 
        name = 'profile'
    else:
        if '.' in out:
            name = out[:out.index(".")]
        else:
            name = out

    fig.print_fig(name)

