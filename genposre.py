#!/usr/bin/env python3
# coding: utf-8

# #### Last update: Nov. 8 2023
# #### It generates position restraints file for Gromacs.

# beer-ware licence
# oanca.gabriel@gmail.com

# For download and updates, vizit or clone:
#     https://github.com/gabrieloanca/gmxtools.git  
#     git@github.com:gabrieloanca/gmxtools.git
# For suggestions, reporting buggs or for any assistance write to oanca.gabriel@gmail.com

import sys, getopt

def show_help():
    print('''
genposre.py [-h] [-q qmatoms.dat] [-i posre.itp] [-r1 QM atoms constraints] [-r2 MM atoms constraints] [-o posre<r1><r2>.itp]

Options:
 -h                                    show this help and exit
 -q or --qmatoms  [qmatoms.dat]        the same file used to give region 1 parameters in gmx4evb.py
 -i or --input    [posre.itp]          posre.itp file generate by pdb2gmx tool of gromacs
 -o or --output   [posre<r1><r2>.itp]  the newly generated constraints file
 -1 or --r1       [0.5]                the constraining force constant in  kcal/mol/Å^2 for atoms in region 1
 -2 or --r2       [0.0]                the constraining force constant in  kcal/mol/Å^2 for atoms in region 2
''')

def get_args():
    arg_list = sys.argv[1:]
    options="hq:i:o:1:2:"
    long_options=["help","qmatoms=","input=","output=","r1=","r2="]
    opts, args= getopt.getopt(arg_list, options, long_options)
    qfile = 'qmatoms.dat'
    inp = 'posre.itp'
    out = ''
    c1 = 0.5
    c2 = 0.0
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            show_help()
            sys.exit()
        elif opt in ('-q', '--qmatoms'):
            qfile = str(arg)
        elif opt in ('-i', '--input'):
            inp = str(arg)
        elif opt in ('-o', '--output'):
            out = str(arg)
        elif opt in ('-1', '--r1'):
            c1 = float(arg)
        elif opt in ('-2', '--r2'):
            c2= float(arg)           
    
    #convert to kJ/mol/nm2
    con1 = c1 / 0.239006 * 100
    con2 = c2 / 0.239006 * 100

    return qfile, inp, out, con1, con2, c1, c2

def qatoms(qfile):
    r1 = [] # list of pdb indeces for region 1 atoms

    with open(qfile) as f:
        data = f.read().split("\n")

    for i,l in enumerate(data):
        if ';' in l:
            ind = l.index(';')
            l = l[:ind]
        elif '#' in l:
            ind = l.index('#')
            l = l[:ind]

        line = l.strip().split()
        if (not line):
            continue
        elif ('[atoms]' in ''.join(line)):
            at = True
        elif ('[bonds]' in ''.join(line)):
            at = False
        elif ('[bcon]' in ''.join(line)):
            at = False
        elif ('[soft-core]' in ''.join(line)):
            at = False
        elif ('[soft-pairs]' in ''.join(line)):
            at = False
        elif ('[angles]' in ''.join(line)):
            at = False
        elif ('[torsions]' in ''.join(line)):
            at = False
        elif ('[impropers]' in ''.join(line)):
            at = False
        else:
            if at:
                try:
                    if (line[6]=='1'):
                        r1.append(int(line[0]))
                except:
                    pass
    return r1

def genposre(inp, r1, con1, con2):
    with open(inp) as f:
        data = f.read().split("\n")

    for line in data:
        l = line.strip().split()
        if "[position_restraints]" in ''.join(l):
            ind = data.index(line)
            break

    for i in range(ind+2, len(data)-1):
        line = data[i]
        line = line.split()
        if int(line[0]) in r1:
            newline = f'{line[0]:>6}     {line[1]} {con1:>12.4f} {con1:12.4f} {con1:12.4f} {con1:>12.4f} {con1:12.4f} {con1:12.4f}'
        else:
            newline = f'{line[0]:>6}     {line[1]} {con2:>12.4f} {con2:12.4f} {con2:12.4f} {con2:>12.4f} {con2:12.4f} {con2:12.4f}'
        data[i] = newline

    return data

def write(data, inp, out, con1, con2):
    if not out:
        k1 = k2 = ''

        if (1 > con2 > 0) and not (con1 == con2):
            con2 = list(str(con2))
            con2.pop(con2.index('.'))
            k2 = ''.join(con2)

        if con1 >= 1:
            k1 = str(int(con1))
        elif 1 > con1 > 0:
            con1 = list(str(con1))
            con1.pop(con1.index('.'))
            k1 = ''.join(con1)

        base = list(inp)
        base = base[:base.index('.')]
        base = ''.join(base)
        out = base + k1 + k2 + '.itp'

    with open(out , "w") as f:
        for l in data:
            f.write(l + '\n')


if __name__ == "__main__":
    try:
        qfile, inp, out, con1, con2, c1, c2 = get_args()
    except:
        if (sys.argv[1:]) and (sys.argv[1:][0] in ('-h', '--help')):
            pass
        else:
            print('Usage: genposre.py [-q qmatoms.dat] [-i posre.itp] [--r1 QM atoms constraints] [--r2 MM atoms constraints] [-o posre<r1><r2>.itp]')
            print('For help type "genposre.py -h"')
        sys.exit()

    try:
        r1 = qatoms(qfile)
    except:
        print(f'{qfile} file could not be read')
        sys.exit()

    try:
        data = genposre(inp, r1, con1, con2)
    except:
        print(f'{inp} file could not be read')
        sys.exit()

    try:
        write(data, inp, out, c1, c1)
    except:
        print("Output file could not be written")
        sys.exit()

