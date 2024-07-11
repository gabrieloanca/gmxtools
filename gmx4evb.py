#!/usr/bin/env python3
# coding: utf-8

# #### Last update: Apr. 8 2024
# #### It builds topologies, one for each FEP frame, for an EVB simulation with GROMACS.

# beer-ware licence
# oanca.gabriel@gmail.com

# For download and updates, vizit or clone:
#     https://github.com/gabrieloanca/gmxtools.git
#     git@github.com:gabrieloanca/gmxtools.git
# For suggestions, reporting bugs or for any assistance write to oanca.gabriel@gmail.com


import sys, os 
import copy
import math as m
import argparse
from datetime import datetime

def get_args():
    parser = argparse.ArgumentParser(epilog='''
gmx4evb.py reads ffld2gmx.py .opls files and generates topologies for an EVB simulation with GROMACS.
RS and PS residue names should match the prefix from the ffld2gmx.py generated files that correspond to
the moieties in RS and PS states. More that one residue can be passed to each state. In case when more
residues with the same name are present in the same state, they must be mentioned that many times.
''')
    parser.add_argument("-f", "--frames", help="number of FEP windows", required=True, type=int)
    parser.add_argument("-q", "--qmatoms", help="QM atoms file constaining atom types, charges, soft-core repulsions\
                        and some other bonding parameters for the QM atoms (default: qmatoms.dat)", required=False, default="qmatoms.dat")
    parser.add_argument("-t", "--topology", help="Gromacs generated topology (default: topol.top)", required=False, default="topol.top")
    parser.add_argument("-r", "--reactants", nargs='*', help="list of reacting moieties", required=True, type=list)
    parser.add_argument("-p", "--products", nargs='*', help="list of product moieties", required=True, type=list)
    parser.add_argument("--cutoff", help="cut-off distance (default: 10 Å)", required=False, type=float, default=1.0)
    parser.add_argument("--precision", help="choose between 'single' or 'double' precision, depending on your GROMACS installation (default: single)", required=False, type=str, default='single')
    args = parser.parse_args()

    rs , ps = [], []
    for i in args.reactants:
        rs.append(''.join(i))
    for i in args.products:
        ps.append(''.join(i))

    return args.frames, args.qmatoms, args.topology, rs, ps, args.cutoff, args.precision


def filelist(state, top=''):
    files = {}
    for i in state:
        files[i] = []
    
    flist = os.listdir()
    for f in flist:
        if '.opls' in f:
            for pre in state:
                if f.startswith(pre):
                    if not (f in files[pre]):
                        files[pre].append(f)    

    return files

def read_qm(qm, rs, ps):
    q1, q2 = [], []         # a list of atom-types
    du = {}                 # holds dummy atom {pdb_index: dummy_type, ...}
    soft_at = {}            # holds soft-core {pdb_index: (A, beta), ...}
    charges = {}            # {at#: (ch_A, ch_B, A(SC), beta(SC)), ...}
    bevb, bconstr, soft_pairs, aevb, torevb, impevb = [], [], [], [], [], []
    at, bnds, bcon, excl, sc, psc, angs, tors, imps = False, False, False, False, False, False, False, False, False
    
    #define atom type dictionaries (for finding pdb indexes)
    def q2pdb(state, qpdb, st):
        qpdb[st] = []
        for i in state:
            if (state.count(i)) > 1:
                qpdb[st].append({})
            elif state.count(i) == 1 and (not (qpdb[st])):
                qpdb[st].append({})
    qpdb = {} # e.g.: {'rs': {'h2o': [o, h1, h2]}}
    q2pdb(rs, qpdb, 'rs')
    q2pdb(ps, qpdb, 'ps')
       
    with open(qm) as f:
        data = f.read().split("\n")

    for i,l in enumerate(data):
        if ';' in l:
            ind = l.index(';')
            l = l[:ind]
        elif '#' in l:
            ind = l.index('#')
            l = l[:ind]

        line = l.strip().split()
        try:
            if (not line):
                continue
            elif ('[atoms]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = True, False, False, False, False, False, False, False
            elif ('[bonds]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = False, True, False, False, False, False, False, False
            elif ('[bcon]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = False, False, True, False, False, False, False, False
            elif ('[soft-core]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = False, False, False, True, False, False, False, False
            elif ('[soft-pairs]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = False, False, False, False, True, False, False, False
            elif ('[angles]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = False, False, False, False, False, True, False, False
            elif ('[torsions]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = False, False, False, False, False, False, True, False
            elif ('[impropers]' in ''.join(line)):
                at, bnds, bcon, sc, psc, angs, tors, imps = False, False, False, False, False, False, False, True
            else:
                if at:
                    # 'ind' is checking if there is more than one residue with the same name
                    ind = 0
                    while (line[1] in qpdb['rs'][ind].keys()):
                        ind +=1
                    qpdb['rs'][ind][line[1]] = (line[0], int(line[-1]))
    
                    ind = 0
                    while (line[3] in qpdb['ps'][ind].keys()):
                        ind += 1
                    # line[6] also needed for ps in atom_list()
                    qpdb['ps'][ind][line[3]] = (line[0], int(line[-1]))
                    
                    if not (line[1] in q1):
                        q1.append(line[1]) # [rs_type]
                    if not (line[3] in q2):
                        q2.append(line[3]) # [ps_type]
                    charges[line[0]] = (float(line[2]), float(line[4]))
        
                    if line[-1] == '1':
                        if len(line) == 7:
                            du[line[0]] = line[5]
                        else:
                            sys.exit()

                elif sc:
                    soft_at[line[0]] = (float(line[1]), float(line[2])) # {'pdb_index': (A, beta)}
                elif bnds:
                    bevb.append(line)
                elif bcon:
                    bconstr.append(l)
                elif psc:
                    soft_pairs.append(line)
                elif angs:
                    aevb.append(line)
                elif tors:
                    torevb.append(line)
                elif imps:
                    impevb.append(line)
        except:
            print(f"\nTere is an error in {qm} file at line {i+1}.")
            sys.exit()

    return qpdb, q1, q2, du, charges, soft_at, bevb, bconstr, soft_pairs, aevb, torevb, impevb

# build a list of [pdb#, atom_type(PS), charge(PS)]
# used only in write_top()
# I use 'du' because I have here all pdb atom numbers
def atom_list(charges, qpdb):
    atoms = {}
    for j in qpdb['rs']:
        for at in j.keys():
            if j[at][1] == 1: 
                ind = j[at][0]
                # I need the charges later, when adding state B to atoms in topology
                qA = charges[ind][0]
                qB = charges[ind][1]
                atoms[ind] = [qA, at, qB]
    for j in qpdb['ps']:
        for at in j.keys():
            if j[at][1] == 1:
                ind = j[at][0]
                atoms[ind].append(at)
    
    return atoms

def bonds_list(rs, ps, rs_files, ps_files, qpdb, q1, q2, bevb):
    harmonic = 1
    
    rsbf = []      # bond files for RS state
    psbf = []      # bond files for PS state
    rs_bonds = []  # bonds from rsbf files
    ps_bonds = []  # bonds from psbf files
    bonds = []     # bonds in RS and PS paired by their pdb numbers
    bonds_x2y = [] # used for coulomb 1-2 and soft-core

    def get_files(state, files, stbf):
        for i in state:
            for j in files[i]:
                if 'bonds' in j:
                    if not (j in stbf):
                        stbf.append(j) # 'stbf': state_bond_file
    
    # 'stbf' is a state_bond_file
    def get_bonds(files, st_qpdb, q, bonds):
        for file in files:
            with open(file) as f:
                data = f.read().strip().split("\n")
            for line in data:
                line = line.split()
                if line[0] in ['#', ";", "!"]:
                    pass
                elif (line[0] in q) and (line[1] in q):
                    for ql in st_qpdb:
                        if ((line[0] in ql.keys()) and (line[1] in ql.keys())) and (int(ql[line[0]][1]) == int(ql[line[1]][1]) == 1):
                            bnd = [ql[line[0]][0], ql[line[1]][0], line[3], line[4]]
                            if bnd in bonds:
                                break
                            else:
                                bonds.append(bnd)
            del(data)
    
    # check if a bond is present only in rs or only in ps
    # here I used 5 whenever > 2 (i.e., more than one bond away)
    def check_bonds(bonds, bonds_x2y, rs_bonds, ps_bonds):
        for l1 in rs_bonds:
            chk = True
            for l2 in ps_bonds:
                if (l1[0] in l2[:2]) and (l1[1] in l2[:2]):
                    bonds.append([l1[0], l1[1], harmonic, l1[2], l1[3], l2[2], l2[3]])
                    chk = False
                    break
            if chk:
                #bonds.append([l1[0], l1[1], morse, '<b0>', '<D>', '<beta>', '<b0>', 0.0, 0.0])
                bonds.append([l1[0], l1[1], harmonic, l1[2], l1[3], 0.0, 0.0])
                bonds_x2y.append((l1[:2], 25))
    
        # check if a bond is present only in ps
        for l2 in ps_bonds:
            chk = True
            for l1 in rs_bonds:
                if (l2[0] in l1[:2]) and (l2[1] in l1[:2]):
                    chk = False
                    break
            if chk:
                #bonds.append([l2[0], l2[1], morse, '<b0>', 0.0, 0.0, '<b0>', '<D>', '<beta>'])
                bonds.append([l2[0], l2[1], harmonic, 0.0, 0.0, l2[2], l2[3]])
                bonds_x2y.append((l2[:2], 52))
    
    def substitute_bevb(bonds, bevb):
        for m in bevb:
            chk = True
            for i, b in enumerate(bonds):
                if (m[0] in b[:2]) and (m[1] in b[:2]):
                    chk = False
                    bonds[i] = m
                    break
            if chk:
                bonds.append(m)
            
    get_files(rs, rs_files, rsbf)
    get_files(ps, ps_files, psbf)            
            
    get_bonds(rsbf, qpdb['rs'], q1, rs_bonds)
    get_bonds(psbf, qpdb['ps'], q2, ps_bonds)
    
    check_bonds(bonds, bonds_x2y, rs_bonds, ps_bonds)
    
    # insert the bonds defined in qmatoms.dat. We usually substitute forming/breaking bonds with morse
    substitute_bevb(bonds, bevb)
    
    return bonds, bonds_x2y

### atoms in forming/breaking angles must all be in region 1
### because they contribute to EVB bonding energy
def angles_list(rs, ps, rs_files, ps_files, qpdb, q1, q2, aevb):
    func = 1
    
    rsaf = []       # anlge files for RS state
    psaf = []       # angle files for PS state
    rs_angles = []  # angles from rsaf files
    ps_angles = []  # angles from psaf files
    angles = []     # angles in RS and PS paired by their pdb numbers
    angles_x2y = [] # used for coulomb 1-2 and soft-core

    def get_files(state, files, staf):
        for i in state:
            for j in files[i]:
                if 'angles' in j:
                    if not (j in staf):
                        staf.append(j) # 'staf': state_angle_file
    
    def get_angles(files, st_qpdb, q, angles):
        for file in files:
            with open(file) as f:
                data = f.read().strip().split("\n")
            for line in data: 
                line = line.split()
                if line[0] in ['#', ";", "!"]:
                    pass
                elif (line[0] in q) and (line[1] in q) and (line[2] in q):
                    for ql in st_qpdb:
                        if ((line[0] in ql.keys()) and (line[1] in ql.keys()) and (line[2] in ql.keys())) and ((int(ql[line[0]][1]) == int(ql[line[1]][1]) == int(ql[line[2]][1]) == 1)):
                            ang = [ql[line[0]][0], ql[line[1]][0], ql[line[2]][0], line[4], line[5]]
                            if ang in angles:
                                break
                            else:
                                angles.append(ang)           
            del(data)
    
    def check_angles(angles, angles_x2y, rs_angles, ps_angles):
        # check if angles only in rs
        # here I used 5 whenever it is not 3 (i.e., 2 bonds away)
        for l1 in rs_angles:
            chk = True
            for l2 in ps_angles:
                if (l1[:3] == l2[:3]) or (l1[:3] == l2[:3][::-1]):
                    angles.append([l1[0], l1[1], l1[2], func, l1[3], l1[4], l2[3], l2[4]])
                    chk = False
                    break
            if chk:
                angles.append([l1[0], l1[1], l1[2], func, l1[3], l1[4], 0.0, 0.0])
                angles_x2y.append((l1[:3], 35))
    
        # check if angles only in ps
        for l2 in ps_angles:
            chk = True
            for l1 in rs_angles:
                if (l2[:3] == l1[:3]) or (l2[:3] == l1[:3][::-1]):
                    chk = False
                    break
            if chk:
                angles.append([l2[0], l2[1], l2[2], func, 0.0, 0.0, l2[3], l2[4]])
                angles_x2y.append((l2[:3], 53))

    def substitute_aevb(angles, aevb):
        for ang in aevb:
            chk = True
            for i, a in enumerate(angles):
                if (ang[0] in a[:3]) and (ang[1] in a[:3]) and (ang[2] in a[:3]):
                    chk = False
                    angles[i] = ang
                    break
            if chk:
                angles.append(ang)
                
    get_files(rs, rs_files, rsaf)
    get_files(ps, ps_files, psaf)
                
    get_angles(rsaf, qpdb['rs'], q1, rs_angles)
    get_angles(psaf, qpdb['ps'], q2, ps_angles)
    
    check_angles(angles, angles_x2y, rs_angles, ps_angles)
    
    # insert angles defined in qmatoms.dat
    substitute_aevb(angles, aevb)
    
    return angles, angles_x2y

def torsions_list(param, rs, ps, rs_files, ps_files, qpdb, q1, q2, torevb):
    rstf = []         # torsion files for RS state
    pstf = []         # torsion files for PS state
    rs_torsions = []  # torsions from rstf files
    ps_torsions = []  # torsions from pstf files
    torsions = []     # torsions in RS and PS paired by their pdb numbers

    def get_files(state, files, sttf, param):
        for i in state:
            for j in files[i]:
                if param in j:
                    if not (j in sttf):
                        sttf.append(j) # 'sttf': state_torsion_file    
    
    def get_tor(files, st_qpdb, q, torsions):
        for file in files:
            with open(file) as f:
                data = f.read().strip().split("\n")
            for line in data:
                line = line.split()
                if line[0] in ['#', ";", "!"]:
                    pass
                elif (line[0] in q) and (line[1] in q) and (line[2] in q) and (line[3] in q):
                    for ql in st_qpdb:
                        if ((line[0] in ql.keys()) and (line[1] in ql.keys()) and line[2] in ql.keys() and line[3] in ql.keys()) and (int(ql[line[0]][1])+int(ql[line[1]][1])+int(ql[line[2]][1])+int(ql[line[3]][1])<=5):
                            tor = [ql[line[0]][0], ql[line[1]][0], ql[line[2]][0], ql[line[3]][0], line[4], line[5:]]
                            if tor in torsions:
                                break
                            else:
                                torsions.append(tor)
            del(data)        
    
    
    def check_tor(torsions, rs_torsions, ps_torsions):
        if rs_torsions:
            func = int(rs_torsions[0][4])
        elif ps_torsions:
            func = int(ps_torsions[0][4])  
        else:
            func = 0
        
        # if proper torsion
        if func == 3:
            for l1 in rs_torsions:
                chk = True
                for l2 in ps_torsions:
                    if (l1[:4] == l2[:4]) or (l1[:4] == l2[:4][::-1]):
                        torsions.append([l1[0], l1[1], l1[2], l1[3], func, l1[5][0], l1[5][1], l1[5][2], l1[5][3], l1[5][4], l1[5][5], l2[5][0], l2[5][1], l2[5][2], l2[5][3], l2[5][4], l2[5][5]])
                        chk = False
                        break
                if chk:
                    torsions.append([l1[0], l1[1], l1[2], l1[3], func, l1[5][0], l1[5][1], l1[5][2], l1[5][3], l1[5][4], l1[5][5], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) 
                        
            for l2 in ps_torsions:
                chk = True
                for l1 in rs_torsions:
                    if (l2[:4] == l1[:4]) or (l2[:4] == l1[:4][::-1]):
                        chk = False
                        break
                if chk:
                    torsions.append([l2[0], l2[1], l2[2], l2[3], func, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, l2[5][0], l2[5][1], l2[5][2], l2[5][3], l2[5][4], l2[5][5]])

        # if improper torsion
        elif func == 2:
            for l1 in rs_torsions:
                chk = True
                for l2 in ps_torsions:
                    if (l1[:4] == l2[:4]) or (l1[:4] == l2[:4][::-1]):
                        torsions.append([l1[0], l1[1], l1[2], l1[3], func, l1[5][0], l1[5][1], l2[5][0], l2[5][1]])
                        chk = False
                        break
                if chk:
                    torsions.append([l1[0], l1[1], l1[2], l1[3], func, l1[5][0], l1[5][1], l1[5][0], 0.0])
                        
            for l2 in ps_torsions:
                chk = True
                for l1 in rs_torsions:
                    if (l2[:4] == l1[:4]) or (l2[:4] == l1[:4][::-1]):
                        chk = False
                        break
                if chk:
                    torsions.append([l2[0], l2[1], l2[2], l2[3], func, l2[5][0], 0.0, l2[5][0], l2[5][1]])
                    #torsions_x2y.append([l2[0], l2[1], l2[2], l2[3], func, l2[4], 0.0, l2[4], l2[5]])
                    
    def substitute_torevb(torsions, torevb):
        for tor in torevb:
            chk = True
            for i, t in enumerate(torsions):
                if (tor[0] in t[:4]) and (tor[1] in t[:4]) and (tor[2] in t[:4]) and (tor[3] in t[:4]):
                    chk = False
                    torsions[i] = tor
                    break
            if chk:
                torsions.append(tor)            
    
    get_files(rs, rs_files, rstf, param)
    get_files(ps, ps_files, pstf, param)
        
    get_tor(rstf, qpdb['rs'], q1, rs_torsions)
    get_tor(pstf, qpdb['ps'], q2, ps_torsions)

    check_tor(torsions, rs_torsions, ps_torsions)
    
    # insert torsions defined in qmatoms.dat
    substitute_torevb(torsions, torevb)
    
    return torsions

### check for bonds/angles/LJ14 complications as in the following transformations (see pictograms)
### it also checks for bonds/angles/14pairs in 3, 4 and 5 atoms rings
#   *--       *--*
#  /\    ->  /  /
# *-*       *--*
###### and ######
#   *--*      *-*
#  /     ->  /\/
# *--*      *-*
def pairs_list(rs, ps, rs_files, ps_files, qpdb, q1, q2):
    # lists of vdw, bonds, and angle files in RS and PS states
    rsvdwf = []    # vdw files in rs
    rsbf = []      # bonds files in rs
    rsaf = []      # angles files in rs
    psvdwf = []
    psbf = []
    psaf = []
    pairs_x2y = [] # this gets returned
    
    # get_files(state, state_file, state_vdw, state_bond, state_angle)
    def get_files(state, files, xvdwf, xbf, xaf):
        for i in state:
            for j in files[i]:
                if 'vdw' in j:
                    if not (j in xvdwf):
                        xvdwf.append(j)
                elif 'bonds' in j:
                    if not (j in xbf):
                        xbf.append(j)
                elif 'angles' in j:
                    if not (j in xaf):
                        xaf.append(j)
    
    def get_vdws(files, st_qpdb, q):
        vdw={}
        for file in files:
            with open(file) as f:
                data = f.read().strip().split("\n")
            for line in data:
                line = line.split()
                if line[0] in ['#', ";", "!"]:
                    pass
                elif line[0] in q:
                    for ql in st_qpdb:
                        if line[0] in ql.keys():
                            if ql[line[0]][0] in vdw.keys():
                                break
                            else:
                                # {pdb#: (sigma, epsilon)}
                                vdw[ql[line[0]][0]] = (line[6], line[7])
            del(data)  # needed in jupyter-notebook
        return vdw
    
    # gets also bonds between q and non-q atoms
    # that's why I don't keep the bonds from bonds_list()
    def get_bonds(files, st_qpdb, q):
        bonds = []
        for file in files:
            with open(file) as f:
                data = f.read().strip().split("\n")
            for line in data:
                line = line.split()
                if line[0] in ['#', ";", "!"]:
                    pass
                elif (line[0] in q) and (line[1] in q):
                    # st_qpdb is a list of dictionaries
                    for ql in st_qpdb:
                        if (line[0] in ql.keys()) and (line[1] in ql.keys()):
                            bond = (ql[line[0]][0], ql[line[1]][0])
                            if bond in bonds:
                                break
                            else:
                                bonds.append(bond)
            del(data)  # needed in jupyter-notebook
        return bonds
    
    def get_angles(files, st_qpdb, q):
        angles = []
        for file in files:
            with open(file) as f:
                data = f.read().strip().split("\n")
            for line in data:
                line = line.split()
                if line[0] in ['#', ";", "!"]:
                    pass
                elif (line[0] in q) and (line[1] in q) and (line[2] in q):
                    for ql in st_qpdb:
                        if (line[0] in ql.keys()) and (line[1] in ql.keys()) and (line[2] in ql.keys()):
                            ang = (ql[line[0]][0], ql[line[1]][0], ql[line[2]][0])
                            if ang in angles:
                                break
                            else:
                                angles.append(ang)
            del(data)  # needed in jupyter-notebook
        return angles
    
    # gives indexes for 1-4 pairs
    def get_indexes(bonds):
        indexes = []  # atom indexes for 1-4 pairs
        for i, l1 in enumerate(bonds[:-2]):
            at1, at2 = l1[:2]
            for j, l2 in enumerate(bonds[i+1:]):
                at3 = None
                if at2 == l2[0]:
                    at3 = l2[1]
                elif at2 == l2[1]:
                    at3 = l2[0]
                elif at1 == l2[0]:
                    at1 = at2
                    at2 = l2[0]
                    at3 = l2[1]
                elif at1 == l2[1]:
                    at1 = at2
                    at2 = l2[1]
                    at3 = l2[0]
                if at3:
                    # starts from 'i+1' because the bonds may not be listed in order
                    for l3 in bonds[i+1:]:
                        at4 = at0 = None
                        if l3 == l2:
                            continue
                        if at3 == l3[0]:
                            at4 = l3[1]
                        elif at3 == l3[1]:
                            at4 = l3[0]
                        elif at1 == l3[0]:
                            at0 = l3[1]
                        elif at1 == l3[1]:
                            at0 = l3[0]
                        if at4:
                            indexes.append((at1, at4))
                        elif at0:
                            indexes.append((at0, at3))
        return indexes

    # remove pairs that form closed loops (or 3 atoms ring, like in epoxy) - see 1st pictogram
    # iterates for every epoxyde (just in case there is more than one epoxyde)
    # NOTE: when you remove an element of a list, it shifts with respect to 'enumerate' 
    def chk_3loop(indexes):
        chk = True
        while chk:
            for i, j in enumerate(indexes):
                if j[0] == j[1]:
                    indexes.pop(i)
                    break
            chk = False        
    
    # remove redundancies
    def redund(indxs):
        nored = []
        for p in indxs:
            if (p in nored) or (p[::-1] in nored):
                pass
            else:
                nored.append(p)

        return nored
    
    # remove pairs that may also form bonds or angles in 4 or 5 atoms rings
    def purify(nored, angles):
        # checks if pair exists in rs_angles or ps_angles      
        # checking for angles already involves checking over bonds
        pure = []
        for p in nored:
            chk = True
            for a in angles:
                if (p[0] in a) and (p[1] in a):
                    chk = False
                    break
            if chk:
                pure.append(p)
                    
        return pure
    
    # keep pairs that belong only to one of the reactant states
    def sift(a, b):
        holder = []
        for p1 in a:
            chk = True
            for p2 in b:
                if (p1[0] in p2) and (p1[1] in p2):
                    chk = False
                    break
            if chk:
                holder.append(p1)
                
        return holder

    # assign labels to pairs that belong only to one state
    def classify(pairs_x2y, pairs, bonds, angles, state):
        if state == "rs":
            labels = (42, 43, 45)
        elif state == "ps":
            labels = (24, 34, 54)
        for p in pairs:
            # here you must also check the pairs, in case there is a double bond breaking that results in a 2 atoms moiety
            chk_b = False # checks if pair is present in bonds
            chk_a = False # checks if pair is present in angles
            
            for b in bonds:
                if (p[0] in b) and (p[1] in b):
                    chk_b = True
                    break
            for a in angles:
                if (a[0] in p) and (a[2] in p):
                    chk_a = True
                    break 
            
            # e.g. 42 = pair (i.e. 1-4) in RS and bond (i.e. 1-2) in PS
            if chk_b:
                pairs_x2y.append((p, labels[0]))
            elif chk_a:
                pairs_x2y.append((p, labels[1]))
            else:
                pairs_x2y.append((p, labels[2]))
    
    # get lists of files containing vdw, bonds, angles and proper dihedrals
    get_files(rs, rs_files, rsvdwf, rsbf, rsaf)
    get_files(ps, ps_files, psvdwf, psbf, psaf)
    
    # saves all vdw in RS and PS
    rs_vdw = get_vdws(rsvdwf, qpdb['rs'], q1)
    ps_vdw = get_vdws(psvdwf, qpdb['ps'], q2)
    
    # don't use bonds and angles from bond_list() or angle_list() because they don't have region 2 atoms
    rs_bonds = get_bonds(rsbf, qpdb['rs'], q1)
    ps_bonds = get_bonds(psbf, qpdb['ps'], q2)
    
    # saves all angles in RS and PS
    rs_angles = get_angles(rsaf, qpdb['rs'], q1)
    ps_angles = get_angles(psaf, qpdb['ps'], q2)
    
    # get indexes for 1-4 pairs in RS and PS
    rs_indexes = get_indexes(rs_bonds)
    ps_indexes = get_indexes(ps_bonds)
    
    # remove pairs from 3 atoms rings
    chk_3loop(rs_indexes)
    chk_3loop(ps_indexes)
    
    # remove redundant pairs
    rs_nored = redund(rs_indexes)
    ps_nored = redund(ps_indexes)
    
    # do this before checking for RS vs. PS
    rs_pure = purify(rs_nored, rs_angles)
    ps_pure = purify(ps_nored, ps_angles)
    
    # get pairs that belong only to one of the states
    rs_only = sift(rs_pure, ps_pure)
    ps_only = sift(ps_pure, rs_pure)
    
    # assign labels to pairs that belong only to one state
    classify(pairs_x2y, rs_only, ps_bonds, ps_angles, "rs")
    classify(pairs_x2y, ps_only, rs_bonds, rs_angles, "ps")
    
    return pairs_x2y, rs_vdw, ps_vdw

# checks if 25 and 52 bonds, and 35 and 53 angles are present in pairs as 24 or 42 and 34 and 43
# it also checks for 23 and 32
def coulomb_list(bonds_x2y, angles_x2y, pairs_x2y, bonds, angles):
    # saves bonds_x2y and angles_x2y in case you'll change the program to write out only the bonds,
    # angles and torsions that form or break (i.e., B state will be written only in case it changes; 
    # in case it stays the same, state B will not be written out.)
    
    # bonds_solo and angle_solo keeps only bonds/angles that break to full NB interaction (1-5) in one of the states
    bonds_solo = copy.deepcopy(bonds_x2y)
    angles_solo = copy.deepcopy(angles_x2y)
    at_ad = []     # list of donor-acceptor atoms
    
    # Delets 24, 34, 42 and 43 from bonds_x2y and angles_x2y
    # which are already given in pairs_x2y
    for p in pairs_x2y:
        for i, b in enumerate(bonds_solo):
            if (p[0][0] in b[0]) and (p[0][1] in b[0]):
                bonds_solo.pop(i)
                break        
        for j, a in enumerate(angles_solo):
            if (a[0][0] in p[0]) and (a[0][2] in p[0]):
                angles_solo.pop(j)
                break 

    # checks if 25 and 52 are not actually 23 and 32
    chk = True
    while chk:
        chk = False
        for i, b in enumerate(bonds_solo):
            for a in angles:
                if (a[0] in b[0]) and (a[2] in b[0]):
                    if b[1] == 25 and a[7]:
                        bonds_solo.pop(i)
                        chk = True
                        break
                    elif b[1] == 52 and a[5]:
                        bonds_solo.pop(i)
                        chk = True
                        break
            if chk:
                break

    # checks if 35 and 53 are not actually 32 and 23
    chk = True
    while chk:
        chk = False
        for i, a in enumerate(angles_solo):
            for b in bonds:
                if (a[0][0] in b[:2]) and (a[0][2] in b[:2]):
                    if b[2] == 3: # if Morse, the list is larger then if it is harmonic
                        if a[1] == 35 and b[7]:    # b[7] checks for k_Morse(PS)
                            angles_solo.pop(i)
                            chk = True
                            break
                        elif a[1] == 53 and b[4]:  # b[4] checks for k_Morse(RS)
                            angles_solo.pop(i)
                            chk = True
                            break
                    elif b[2] == 1: # if harmonic
                        if a[1] == 35 and b[6]:    # b[6] checks for k_harmonic(PS)
                            angles_solo.pop(i)
                            chk = True
                            break
                        elif a[1] == 53 and b[4]:  # b[4] checks for k_harmonic(RS)
                            angles_solo.pop(i)
                            chk = True
                            break
                        
            if chk:
                break

    # returns a list with donor-acceptor atoms
    # it takes as argument bonds_x2y, not bonds_solo from which 2-4 pairs are canceled (they're present in pairs_x2y)
    # and if a soft-core is given, it won't delete its 1-4 vdW
    def donor_acceptor(bonds_x2y):
        for i, b1 in enumerate(bonds_x2y[:-1]):
            for b2 in bonds_x2y[i+1:]:
                if (b1[1] == 25 and b2[1] == 52) or (b1[1] == 52 and b2[1] == 25):
                    if (b1[0][0] in b2[0]) or b1[0][1] in b2[0]:
                        if (b1[0][0] == b2[0][0]):
                            at_ad.append([b1[0][1], b2[0][1]])
                        elif (b1[0][0] == b2[0][1]):
                            at_ad.append([b1[0][1], b2[0][0]])
                        elif (b1[0][1] == b2[0][0]):
                            at_ad.append([b1[0][0], b2[0][1]])
                        elif (b1[0][1] == b2[0][1]):
                            at_ad.append([b1[0][0], b2[0][0]])
        return at_ad

    at_ad = donor_acceptor(bonds_x2y)

    return bonds_solo, angles_solo, at_ad

def soft_core(at_ad, bonds_solo, angles_solo, soft_at, soft_pairs):
    func = 9
    soft = []
    betas = []
    init_pairs = []
    
    try:
        for p in soft_pairs:
            init_pairs.append((p[0],p[1]))
            betaA, betaB = float(p[3]), float(p[5])

            if betaA == betaB:
                if not (betaA in betas):
                    betas.append(betaA)
                soft.append([p[0], p[1], func, betas.index(betaA), float(p[4]), betas.index(betaA), float(p[6]), betaA])
            else:
                if float(p[4]):
                    if not (betaA in betas):
                        betas.append(betaA)
                    soft.append([p[0], p[1], func, betas.index(betaA), float(p[4]), betas.index(betaA), float(p[6]), betaA])
                if float(p[6]):
                    if not (betaB in betas):
                        betas.append(betaB)
                    soft.append([p[0], p[1], func, betas.index(betaB), float(p[4]), betas.index(betaB), float(p[6]), betaB])     
    except:
        pass
    
    for b in bonds_solo:
        at1, at2 = b[0]
        if not (((at1,at2) in init_pairs) or ((at2,at1) in init_pairs)):
            if (at1 in soft_at.keys()) and (at2 in soft_at.keys()):
                A = soft_at[at1][0]*soft_at[at2][0]
                beta = (soft_at[at1][1]*soft_at[at2][1])**(1/2)
                if not (beta in betas):
                    betas.append(beta)   
                i = betas.index(beta)
        
                if (b[1] == 24) or (b[1] == 25):
                    soft.append([at1, at2, func, i, 0.0, i, A, beta])
                elif (b[1] == 52) or (b[1] == 42):
                    soft.append([at1, at2, func, i, A, i, 0.0, beta])
            
    for a in angles_solo:
        at1, at2, at3 = a[0]
        if not (((at1, at3) in init_pairs) or ((at3, at1) in init_pairs)):
            if (at1 in soft_at.keys()) and (at3 in soft_at.keys()):
                A = soft_at[at1][0]*soft_at[at3][0]
                beta = (soft_at[at1][1]*soft_at[at3][1])**(1/2)
                if not (beta in betas):
                    betas.append(beta)
                i = betas.index(beta)
                
                if (a[1] == 35) or (a[1] == 34):
                    soft.append([at1, at3, func, i, 0.0, i, A, beta])
                elif (a[1] == 53) or (a[1] == 43):
                    soft.append([at1, at3, func, i, A, i, 0.0, beta])

    for p in at_ad:
        p1, p2 = p
        if not (((p1, p2) in init_pairs) or ((p2, p1) in init_pairs)):
            if (p1 in soft_at.keys()) and (p2 in soft_at.keys()):
                A = soft_at[p1][0]*soft_at[p2][0]
                beta = (soft_at[p1][1]*soft_at[p2][1])**(1/2)
                if not (beta in betas):
                    betas.append(beta)
                i = betas.index(beta)
                
                soft.append([p1, p2, func, i, A, i, A, beta])  # donor-acceptor

    return soft, betas

# writes tabulated tables for soft-core potential
def gen_tables(betas, cutoff, prec):
    cutoff += 1 # extends the tables by 1 more nm
    delr = 0.002
    if prec == 'double':
        delr = 0.0005

    for i, beta in enumerate(betas):
        r = 0
        file = open(f'topologies/table_b{i}.xvg', "w")
        while r <= cutoff:
            file.write(f"{r:0<10.4f}       {m.exp(-beta*r):.9E}       {beta*m.exp(-beta*r):.9E}\n")
            r += delr
        file.close()

# gives 1-4 and 1-5 non bonding interactions that are excluded because of the
# forming bonds of type 3 (Morse) in PS state which generates exclusions
# The interaction between these pairs of atoms will be passed under [pairs] and [pairs_nb]
def restore_nb(bonds_solo, angles_solo):
    # nb55 holds the indices for atoms two or three bonds away over the forming-breaking bonds;
    # which gets excluded by grompp because of the (morse) bond with k=0.0 in one of the states.
    nb55 = []
    # get 1-3 pairs
    for i, b1 in enumerate(bonds_solo[:-1]):
        at1, at2 = b1[0]
        for b2 in bonds_solo[i+1:]:
            if not (b1[1] == b2[1]):
                at3, at4 = b2[0]
                if at1 == at3:
                    nb55.append((at2, at4))
                elif at1 == at4:
                    nb55.append((at2, at3))
                elif at2 == at3:
                    nb55.append((at1, at4))
                elif at2 == at4:
                    nb55.append((at1, at3))

    for bond in bonds_solo:
        b1, b2 = bond[0]
        for angle in angles_solo:
            a1, a3 = angle[0][0], angle[0][2]
            if ((bond[1] == 25) and (angle[1] == 53)):
                if (b1 == a1):
                    nb55.append((b2, a3))
                elif (b2 == a1):
                    nb55.append((b1, a3))
                elif (b1 == a3):
                    nb55.append((b2, a1))
                elif (b2 == a3):
                    nb55.append((b1, a1))
            if ((bond[1] == 52) and (angle[1] == 35)):
                if (b1 == a1):
                    nb55.append((b2, a3))
                elif (b2 == a1):
                    nb55.append((b1, a3))
                elif (b1 == a3):
                    nb55.append((b2, a1))
                elif (b2 == a3):
                    nb55.append((b1, a1))
            
    return nb55

# this method is invoked for every lambda, when generating the topology
# it gives the NB pairs to be filled in [pairs_nb]
# C6 = 4*epsilon*sigma**6; C12 = 4*epsilon*sigma**12 (GROMACS manual-2020.3, pg 394)
# in [pairs_nb] sigma does not change, epsilon=(1-l)*RS+l*PS
def pairs_nb(at_ad, soft, nb55, bonds_solo, angles_solo, pairs_x2y, charges, rs_vdw, ps_vdw, l): # l = lambda = 0..1
    func = 1
    nb_list = []   # for [pairs_nb]
    ps_list = []   # for [pairs]
    # nb params are like {at#: (ch_A, ch_B, A(SC), beta(SC))}
    
    # add Coulomb. No vdW; it has soft core instead !!!
    for b in bonds_solo:
        at1, at2 = b[0]
        # NB in state B
        if b[1] == 25:
            q1B = float(charges[at1][1])
            q2B = float(charges[at2][1])
            qtotB = q1B*q2B
            
            qAB = qtotB*l
            nb_list.append([at1, at2, func, qAB, 1, '  ;2 -> 5'])
            
        # NB in state A
        elif b[1] == 52:
            q1A = float(charges[at1][0])
            q2A = float(charges[at2][0])
            qtotA = q1A*q2A
            
            qAB = qtotA*(1-l)
            nb_list.append([at1, at2, func, qAB, 1, '  ;5 -> 2'])
    
    # add Coulomb. No vdW; it has soft core instead !!!       
    for a in angles_solo:
        at1, at2, at3 = a[0]
        # NB in state B
        if a[1] == 35:
            q1B = float(charges[at1][1])
            q3B = float(charges[at3][1])
            qtotB = q1B*q3B
            
            qAB = qtotB*l
            nb_list.append([at1, at3, func, qAB, 1, '  ;3 -> 5'])            
        
        # NB in state A
        elif a[1] == 53:
            q1A = float(charges[at1][0])
            q3A = float(charges[at3][0])
            qtotA = q1A*q3A            
            
            qAB = qtotA*(1-l)
            nb_list.append([at1, at3, func, qAB, 1, '  ;5 -> 3'])            
            
    for p in pairs_x2y:
        at1, at4 = p[0]
        tick = True
        if ([at1, at4] in at_ad) or ([at4, at1] in at_ad):
            for s in soft:
                if ([at1, at4] == s[:2]) or ([at4, at1] == s[:2]):
                    tick = False
                    break
                    
            if tick:
                if l == 0 :
                    print(f"The donor acceptor atoms {at1} and {at4} does not have soft-core repulsion.")
                    print("Since they are 4 bonds away, 1-4 van der Waals have been added instead.")
                    print("If you want to avoid adding 1-4 van der Waals, then add a null soft-core.")
                    print()
            else:
                # 0.5*NB in state B
                if (p[1] == 24) or (p[1] == 34):
                    #q1A = float(charges[at1][0])
                    #q4A = float(charges[at4][0])
                    #qtotA = q1A*q4A            
                    
                    q1B = float(charges[at1][1])
                    q4B = float(charges[at4][1])
                    qtotB = q1B*q4B
                    
                    qAB = 0.5*qtotB*l
                    nb_list.append([at1, at4, func, qAB, 1, '  ;2/3 -> 4; soft core']) 
                
                # 0.5*NB in state A
                elif (p[1] == 42) or (p[1] == 43):
                    q1A = float(charges[at1][0])
                    q4A = float(charges[at4][0])
                    qtotA = q1A*q4A            
                    
                    #q1B = float(charges[at1][1])
                    #q4B = float(charges[at4][1])
                    #qtotB = q1B*q4Bs
                    
                    qAB = 0.5*qtotA*(1-l)
                    nb_list.append([at1, at4, func, qAB, 1, '  ;4 -> 2/3; soft core'])            
                    
                #0.5*NB in state A + 1*NB in state B
                elif p[1] == 45:
                    q1A = float(charges[at1][0])
                    q4A = float(charges[at4][0])
                    qtotA = q1A*q4A            
                    
                    q1B = float(charges[at1][1])
                    q4B = float(charges[at4][1])
                    qtotB = q1B*q4B
                    
                    qAB = 0.5*qtotA*(1-l) + qtotB*l
                    nb_list.append([at1, at4, func, qAB, 1, '  ;4 -> 5; soft core'])            
                
                # full NB in state A + 0.5 * NB in state B
                elif p[1] == 54:
                    q1A = float(charges[at1][0])
                    q4A = float(charges[at4][0])
                    qtotA = q1A*q4A            
                    
                    q1B = float(charges[at1][1])
                    q4B = float(charges[at4][1])
                    qtotB = q1B*q4B
                    
                    qAB = qtotA*(1-l) + 0.5*qtotB*l
                    nb_list.append([at1, at4, func, qAB, 1, '  ;5 -> 4; soft core']) 
                
        #else:
        if tick:
            # 0.5*NB in state B
            if (p[1] == 24) or (p[1] == 34):
                q1A = float(charges[at1][0])
                q4A = float(charges[at4][0])
                qtotA = q1A*q4A            
                
                q1B = float(charges[at1][1])
                q4B = float(charges[at4][1])
                qtotB = q1B*q4B
                
                #s1A = float(rs_vdw[at1][0])
                #s4A = float(rs_vdw[at4][0])
                #stotA = (s1A*s4A)**0.5
                
                #e1A = float(rs_vdw[at1][1])
                #e4A = float(rs_vdw[at4][1])
                #etotA = (e1A*e4A)**0.5
                
                s1B = float(ps_vdw[at1][0])
                s4B = float(ps_vdw[at4][0])
                stotB = (s1B*s4B)**0.5
                
                e1B = float(ps_vdw[at1][1])
                e4B = float(ps_vdw[at4][1])
                etotB = 0.5 * (e1B*e4B)**0.5
                
                qAB = - 0.5*qtotA*(1-l)
                ps_list.append([at1, at4, func, 1.0, 0.0, stotB, etotB, '  ;2/3 -> 4'])
                nb_list.append([at1, at4, func, qAB, 1, '  ;2/3 -> 4']) 
            
            # 0.5*NB in state A
            elif (p[1] == 42) or (p[1] == 43):
                q1A = float(charges[at1][0])
                q4A = float(charges[at4][0])
                qtotA = q1A*q4A            
                
                q1B = float(charges[at1][1])
                q4B = float(charges[at4][1])
                qtotB = q1B*q4B
                
                s1A = float(rs_vdw[at1][0])
                s4A = float(rs_vdw[at4][0])
                stotA = (s1A*s4A)**0.5
                
                e1A = float(rs_vdw[at1][1])
                e4A = float(rs_vdw[at4][1])
                etotA = 0.5 * (e1A*e4A)**0.5
                
                #s1B = float(ps_vdw[at1][0])
                #s4B = float(ps_vdw[at4][0])
                #stotB = (s1B*s4B)**0.5
                
                #e1B = float(ps_vdw[at1][1])
                #e4B = float(ps_vdw[at4][1])
                #etotB = (e1B*e4B)**0.5
                
                qAB = - 0.5*qtotB*l
                ps_list.append([at1, at4, func, stotA, etotA, 1.0, 0.0, '  ;4 -> 2/3'])
                nb_list.append([at1, at4, func, qAB, 1, '  ;4 -> 2/3'])            
                
            #0.5*NB in state A + 1*NB in state B
            elif p[1] == 45:
                #q1A = float(charges[at1][0])
                #q4A = float(charges[at4][0])
                #qtotA = q1A*q4A            
                
                q1B = float(charges[at1][1])
                q4B = float(charges[at4][1])
                qtotB = q1B*q4B
                
                s1A = float(rs_vdw[at1][0])
                s4A = float(rs_vdw[at4][0])
                stotA = (s1A*s4A)**0.5
                
                e1A = float(rs_vdw[at1][1])
                e4A = float(rs_vdw[at4][1])
                etotA = 0.5 * (e1A*e4A)**0.5
                
                s1B = float(ps_vdw[at1][0])
                s4B = float(ps_vdw[at4][0])
                stotB = (s1B*s4B)**0.5
                
                e1B = float(ps_vdw[at1][1])
                e4B = float(ps_vdw[at4][1])
                etotB = (e1B*e4B)**0.5
                
                qAB = 0.5*qtotB*l
                ps_list.append([at1, at4, func, stotA, etotA, stotB, etotB, '  ;4 -> 5'])
                nb_list.append([at1, at4, func, qAB, 1, '  ;4 -> 5'])            
            
            # full NB in state A + 0.5 * NB in state B
            elif p[1] == 54:
                q1A = float(charges[at1][0])
                q4A = float(charges[at4][0])
                qtotA = q1A*q4A            
                
                #q1B = float(charges[at1][1])
                #q4B = float(charges[at4][1])
                #qtotB = q1B*q4B
                
                s1A = float(rs_vdw[at1][0])
                s4A = float(rs_vdw[at4][0])
                stotA = (s1A*s4A)**0.5
                
                e1A = float(rs_vdw[at1][1])
                e4A = float(rs_vdw[at4][1])
                etotA = (e1A*e4A)**0.5
                
                s1B = float(ps_vdw[at1][0])
                s4B = float(ps_vdw[at4][0])
                stotB = (s1B*s4B)**0.5
                
                e1B = float(ps_vdw[at1][1])
                e4B = float(ps_vdw[at4][1])
                etotB = 0.5 * (e1B*e4B)**0.5
                
                qAB = 0.5*qtotA*(1-l)
                ps_list.append([at1, at4, func, stotA, etotA, stotB, etotB, '  ;5 -> 4'])
                nb_list.append([at1, at4, func, qAB, 1, '  ;5 -> 4'])            
            
    # restore the excluded pairs due to forming/breaking bonds which generate exclusions
    for p in nb55:
        at1, at2 = p
        tick = True
        if ([at1, at2] in at_ad) or ([at2, at1] in at_ad):
            for s in soft:
                if ([at1, at2] == s[:2]) or ([at2, at1] == s[:2]):
                    tick = False
                    break
                    
            if tick:
                if l == 0:
                    print(f"The donor acceptor atoms {at1} and {at2} does not have soft-core repulsion and full der Waals was added instead.")
                    print("If you want to avoid the van der Waals interaction, then add a null soft-core (i.e., pre-exponential = 0).")
                    print()
            else:
                q1A = float(charges[at1][0])
                q2A = float(charges[at2][0])
                qtotA = q1A*q2A            
                
                q1B = float(charges[at1][1])
                q2B = float(charges[at2][1])
                qtotB = q1B*q2B
                
                qAB = qtotA*(1-l) + qtotB*l
                nb_list.append([at1, at2, func, qAB, 1, '  ;5 -> 5; soft core'])                
        
        if tick:
            q1A = float(charges[at1][0])
            q2A = float(charges[at2][0])
            qtotA = q1A*q2A            
            
            q1B = float(charges[at1][1])
            q2B = float(charges[at2][1])
            qtotB = q1B*q2B
            
            s1A = float(rs_vdw[at1][0])
            s2A = float(rs_vdw[at2][0])
            stotA = (s1A*s2A)**0.5
            
            e1A = float(rs_vdw[at1][1])
            e2A = float(rs_vdw[at2][1])
            etotA = (e1A*e2A)**0.5
            
            s1B = float(ps_vdw[at1][0])
            s2B = float(ps_vdw[at2][0])
            stotB = (s1B*s2B)**0.5
            
            e1B = float(ps_vdw[at1][1])
            e2B = float(ps_vdw[at2][1])
            etotB = (e1B*e2B)**0.5
            
            qAB = 0.5*qtotA*(1-l) + 0.5*qtotB*l
            ps_list.append([at1, at2, func, stotA, etotA, stotB, etotB, '  ;5 -> 5'])
            nb_list.append([at1, at2, func, qAB, 1, '  ;5 -> 5'])
        
    return nb_list, ps_list

### remove pairs_4x from topology - they will be inserted in [pairs_nb]
### !!! LJ14 chages from A to B when you change the atom_types
def write_top(at_ad,top,atoms,bonds,angles,torsions,impropers,soft,bonds_solo,angles_solo,pairs_x2y,charges,rs_vdw,ps_vdw,feps,bconstr,nb55):
    # define here some parameters only once
    user = os.environ.get('USER')
    date = datetime.now()
    date = str(date.year)+'-'+str(date.month)+'-'+str(date.day)+' '+str(date.hour)+':'+str(date.minute)+':'+str(date.second)
    with open(top) as f:
        data = f.read().split("\n")
        
    #### holds the index of the inserted lines from bonds, pairs, anlges, torsions, and impropers.
    ###b, p, a, t, it = [], [], [], [], []
    # trigers for atoms, bonds, pairs, angles, torsions
    # impropers are treated together with torsions
    ats, bnds, prs, ang, tor, smth = False, False, False, False, False, False
    

    ### 'smth' stays for 'someting else'
    alist = list(atoms.keys())
    for i, line in enumerate(data):
        if "; Include Position restraint file" in line:
            start = i - 1
            break
        if (not line.strip()) or (line.strip()[0] == ';'): # for empty or commented lines
            continue
        elif '[atoms]' in ''.join(line.strip().split()):
            ats, bnds, prs, ang, tor, smth = True, False, False, False, False, False
            continue
        elif '[bonds]' in ''.join(line.strip().split()):
            ats, bnds, prs, ang, tor, smth = False, True, False, False, False, False
            continue
        elif '[pairs]' in ''.join(line.strip().split()):
            ats, bnds, prs, ang, tor, smth = False, False, True, False, False, False
            continue
        elif '[angles]' in ''.join(line.strip().split()):
            ats, bnds, prs, ang, tor, smth = False, False, False, True, False, False
            continue
        elif '[dihedrals]' in ''.join(line.strip().split()):
            ats, bnds, prs, ang, tor, smth = False, False, False, False, True, False
            continue
        elif ('[' in ''.join(line.strip().split())) and (']' in ''.join(line.strip().split())):
            ats, bnds, prs, ang, tor, smth = False, False, False, False, False, True
        else:    
            if ats:
                l = line.split()
                # atlist is define before the this loop
                if l[0] in alist:
                    # add [type, chargeB, mass]
                    tpA, chA = atoms[l[0]][1], atoms[l[0]][0]
                    tpB, chB = atoms[l[0]][3], atoms[l[0]][2]
                    # takes care of comments at the end of line
                    if ';' in line:
                        com = line.index(";")
                        tail = line[com:]
                    else:
                        tail = ""
                    line = f' {l[0]:>6} {tpA:>9} {l[2]:>6} {l[3]:>7} {l[4]:>5} {l[5]:>6} {chA:10.6f} {l[7]:>10} {tpB:>8} {chB:10.6f} {float(l[7]):>10}   ' + tail
                    data[i] = line
                        
            elif bnds:
                l = line.split()
                for lx in bonds:
                    if (l[0] in lx[:2]) and (l[1] in lx[:2]):
                        newline = '; ' + line
                        #newline = f' {lx[0]:>5} {lx[1]:>5}    {lx[2]}     {lx[3]}     {lx[4]}     {lx[5]}     {lx[6]}'
                        data[i] = newline
                        break
                
            elif prs:
                l = line.split()
                for lx in pairs_x2y:
                    if 40 < lx[1] < 50:
                        if (l[0] in lx[0]) and (l[1] in lx[0]):
                            newline = '; ' + line
                            #newline = f' {lx[0]:>5} {lx[1]:>5}    {lx[2]}     {lx[3]:.6f}     {lx[4]:.6f}     {lx[5]:.6f}     {lx[6]:.6f}'
                            data[i] = newline
                            break
                            
            elif ang:
                l = line.split()
                for lx in angles:
                    if (l[0] in lx[:3]) and (l[1] in lx[:3]) and (l[2] in lx[:3]):
                        newline = '; ' + line
                        #newline = f' {lx[0]:>5} {lx[1]:>5} {lx[2]:>5}    {lx[3]}     {lx[4]}     {lx[5]}     {lx[6]}     {lx[7]}'
                        data[i] = newline
                        break
            
            elif tor:
                l = line.split()
                for lx in torsions:
                    if (l[0] in lx[:4]) and (l[1] in lx[:4]) and (l[2] in lx[:4]) and (l[3] in lx[:4]):
                        newline = '; ' + line
                        #newline = f' {lx[0]:>5} {lx[1]:>5} {lx[2]:>5} {lx[3]:>5}  {lx[4]}  {lx[5]:>9}  {lx[6]:>9}  {lx[7]:>9}  {lx[8]:>9}  {lx[9]:>9}  {lx[10]:>9}  {lx[11]:>9}  {lx[12]:>9}  {lx[13]:>9}  {lx[14]:>9}  {lx[15]:>9}  {lx[16]:>9}'
                        data[i] = newline
                        break
                        
                for lx in impropers:
                    if (l[0] in lx[:4]) and (l[1] in lx[:4]) and (l[2] in lx[:4]) and (l[3] in lx[:4]):
                        newline = '; ' + line
                        #newline = f' {lx[0]:>5} {lx[1]:>5} {lx[2]:>5} {lx[3]:>5}  {lx[4]}  {lx[5]:>9}  {lx[6]:>9}  {lx[7]:>9}  {lx[8]:>9}  {lx[9]:>9}  {lx[10]:>9}  {lx[11]:>9}  {lx[12]:>9}  {lx[13]:>9}  {lx[14]:>9}  {lx[15]:>9}  {lx[16]:>9}'
                        data[i] = newline
                        break

            elif smth:
                data[i] = line

    
    data.insert(start, '')
    start += 1
    data.insert(start, ';----------------------------------------')
    start += 1
    data.insert(start, '; This section is dedicated to EVB atoms')
    start += 1
    data.insert(start, ';----------------------------------------')
    start += 1
    data.insert(start, '')
    start += 1    
    
    # insert bonds
    chkb = True #check if there are evb bonds
    if bonds:
        chkb = False
        data.insert(start, '[ bonds ]')
        start += 1
        # add harmonic and Morse bonds
        data.insert(start, '; harmonic and Morse bonds')
        start += 1
        for lx in bonds:
            if int(lx[2]) == 3:
                newline = f' {lx[0]:>5} {lx[1]:>5}    {lx[2]}     {lx[3]:>9}     {lx[4]:>9}     {lx[5]:>9}     {lx[6]:>9}     {lx[7]:>9}     {lx[8]:>9}'
                data.insert(start, newline)
                start += 1
            elif (lx[2]) == 1:
                newline = f' {lx[0]:>5} {lx[1]:>5}    {lx[2]}     {lx[3]:>9}     {lx[4]:>9}     {lx[5]:>9}     {lx[6]:>9}'
                data.insert(start, newline)
                start += 1
            
    # add soft-core as tabulated bonds of type 9
    if soft:
        if chkb:
            chkb = False
            data.insert(start, '[ bonds ]')
            start += 1

        data.insert(start, '; soft-core potential')
        start += 1
        for lx in soft:
            newline = f' {lx[0]:>5} {lx[1]:>5}    {lx[2]}    {lx[3]}  {lx[4]:12.2f}    {lx[5]}  {lx[6]:12.2f}  ; beta = {lx[7]:.2f}'
            data.insert(start, newline)
            start += 1
    
    # add bond constraints
    if bconstr:
        if chkb:
            chkb = False
            data.insert(start, '[ bonds ]')
            start += 1

        data.insert(start, '; constraints')
        start += 1
        for lx in bconstr:
            data.insert(start, lx)
            start += 1
    data.insert(start, '')
    start += 1
    
    # insert angles
    if angles:
        data.insert(start, '[ angles ]')
        start += 1
        for lx in angles:
            newline = f' {lx[0]:>5} {lx[1]:>5} {lx[2]:>5}    {lx[3]}     {lx[4]:>9}     {lx[5]:>9}     {lx[6]:>9}     {lx[7]:>9}'
            data.insert(start, newline)
            start += 1
        data.insert(start, '')
        start += 1
    
    # insert torsions
    if torsions:
        data.insert(start, '[ dihedrals ]')
        start += 1
        # add proper dihedrals
        data.insert(start, '; proper dihedrals')
        start += 1
        for lx in torsions:
            newline = f' {lx[0]:>5} {lx[1]:>5} {lx[2]:>5} {lx[3]:>5}  {lx[4]}  {lx[5]:>9}  {lx[6]:>9}  {lx[7]:>9}  {lx[8]:>9}  {lx[9]:>9}  {lx[10]:>9}  {lx[11]:>9}  {lx[12]:>9}  {lx[13]:>9}  {lx[14]:>9}  {lx[15]:>9}  {lx[16]:>9}'
            data.insert(start, newline)
            start += 1

    # add improper dihedrals
    if impropers:
        if not torsions:
            data.insert(start, '[ dihedrals ]')
            start += 1
        data.insert(start, '; improper dihedrals')
        start += 1
        for lx in impropers:
            newline = f' {lx[0]:>5} {lx[1]:>5} {lx[2]:>5} {lx[3]:>5}  {lx[4]}  {lx[5]:>9}  {lx[6]:>9}  {lx[7]:>9}  {lx[8]:>9}'
            data.insert(start, newline)
            start += 1
        data.insert(start, '')
        start += 1
            
    # inserts [exclusions]
    nb_list, ps_list = pairs_nb(at_ad, soft, nb55, bonds_solo, angles_solo, pairs_x2y, charges, rs_vdw, ps_vdw, 1)
    if nb_list:
        data.insert(start, '[ exclusions ]')
        start += 1 
        for lx in nb_list:
            data.insert(start, f' {lx[0]}   {lx[1]}')
            start += 1
        data.insert(start, "")
        start += 1
    
        # insert pairs
        data.insert(start, '[ pairs ]')
        start += 1
        for lx in ps_list:
            newline = f' {lx[0]}   {lx[1]}    {lx[2]}   {lx[3]:>10.6f}  {lx[4]:>10.6f}  {lx[5]:>10.6f}  {lx[6]:>10.6f}'
            try:
                if str(lx[7]):
                    newline = newline+lx[7]
            except:
                pass
            data.insert(start, newline)
            start += 1
        data.insert(start, "")
        start += 1    
    
    # do it here because a top with different [ pairs_nb ] will be built for each FEP frame
    try:
        os.mkdir('topologies')
    except FileExistsError:
        pass
    
    # insert pairs_nb
    data.insert(start, '[ pairs_nb ]')
    start += 1
    nb_start = start
    for i in range(feps):
        nb, _ = pairs_nb(at_ad, soft, nb55, bonds_solo, angles_solo, pairs_x2y, charges, rs_vdw, ps_vdw, i/(feps-1))
        
        ## I use the if/else for not building the topology from the beginning each time
        if i == 0:
            for lx in nb:
                newline = f' {lx[0]:>5} {lx[1]:>5}  {lx[2]}  {lx[3]:>10.6f}  {lx[4]:4.2f}   1.00   0.00'
                try:
                    if str(lx[5]):
                        newline = newline + lx[5]
                except:
                    pass
                data.insert(start, newline)
                start += 1
            #data.insert(start, '')
            #start += 1
        
        elif i > 0:
            for j, lx in enumerate(nb):
                newline = f' {lx[0]:>5} {lx[1]:>5}  {lx[2]}  {lx[3]:>10.6f}  {lx[4]:4.2f}   1.00   0.00'
                try:
                    if str(lx[5]):
                        newline = newline + lx[5]
                except:
                    pass
                data[nb_start+j] = newline
                
        with open(f'topologies/topol_{i:0>3}.top', "w") as f:
            f.write(f'''; Topology for EVB simulation in Gromacs, generated with gmx4evb.py
; User: {user}
; Date: {date}
; For download and updates, vizit or clone:
;     https://github.com/gabrieloanca/gmxtools
;     git@github.com:gabrieloanca/gmxtools
; For suggestions, reporting buggs or for any assistance write to oanca.gabriel@gmail.com
; ---------------------------------------------------------------------------------------
''')
            for l in data:
                f.write(str(l)+'\n')
                
### this function reads topol_000.top and writes out a topology in which all atoms in region 1 will have dummy types
### exclusions will be preserved, soft-core, pair_nb and constraints will be commented out, and bonds, angles,
### torsions and impropers will be set to 0.
### The dummy types must preserve the original bonding types (in ffnonbonded.itp) - this way, the connectivities
### between region 1 and region 2 will be conserved. These dummy types needs to be added to ffnonbonded.itp 
### and atomtypes.atp files of the force field.
def evb_less(top, name, du):  
    ati = list(du.keys())   # holds QM-atom indexes

    with open(f'topologies/{top}') as f:
        data = f.read().split("\n")
    
    new_top = []
    trigger = None # it shows which directive is currently reading
    evb = False    # it shows if it's reading the EVB section
    
    for i, line in enumerate(data):
        if "This section is dedicated to EVB atoms" in line:
            evb = True
            new_top.append(line)
            continue
            
        l = line.strip().split()
        
        # go through triggers
        if ("[atoms]" in ''.join(l)) and (not evb):
            trigger = 'atoms'
            new_top.append(line)
            continue
            
        #elif ("[exclusions]" in ''.join(l)) and (evb):
        #    trigger = 'exclusions'
        #    new_top.append(line)
        #    continue
    
        elif ("[bonds]" in ''.join(l)) and (evb):
            trigger = 'bonds'
            new_top.append(line)
            continue
            
        elif ("soft-core" in line) and (evb):
            trigger = 'soft-core'
            new_top.append(line)
            continue
            
        elif ("constraints" in line) and (evb):
            trigger = 'constraints'
            new_top.append(line)
            continue
            
        elif ("[angles]" in ''.join(l)) and (evb):
            trigger = 'angles'
            new_top.append(line)
            continue
            
        elif ("[dihedrals]" in ''.join(l)) and (evb):
            trigger = 'dihedrals'
            new_top.append(line)
            continue
            
        elif ("[pairs]" in ''.join(l)) and (evb):
            trigger = "pairs"
            new_top.append(line)
            continue

        elif ("[pairs_nb]" in ''.join(l)) and (evb):
            trigger = 'pairs_nb'
            new_top.append(line)
            continue
            
        # This condition must be placed before checking for commented lines
        elif not l:
            trigger = None
            new_top.append(line)
            continue 
            
        # This condition must be placed after all triggers, to get 'soft-core' and 'constraints'
        elif (line.strip()[0] == ';'):
            new_top.append(line)
            continue          
            
        if trigger == "atoms":
            if l[0] in ati:
                dummyA = du[l[0]]
                newline = f'{l[0]:>6} {dummyA:>10} {l[2]:>6} {l[3]:>6} {l[4]:>6} {l[5]:>6}        0.0   {l[7]:>8}'
                new_top.append(newline)
            else:
                new_top.append(line)  
            continue
                
        elif trigger == 'bonds':
            if int(l[2]) == 1:
                newline = f'{l[0]:>6} {l[1]:>5}   1   {l[3]:>8}   0.0'
                new_top.append(newline)
            elif int(l[2]) == 3:
                newline = f'{l[0]:>6} {l[1]:>5}   3   {l[3]:>8}   0.0    {l[5]:>5}'
                new_top.append(newline)
            continue                
        
        elif trigger == 'soft-core':
            newline = ';' + line
            new_top.append(newline)
            continue
            
        elif trigger == 'constraints':
            newline = ";" + line
            new_top.append(newline)
            
        elif trigger == 'angles':    
            newline = f'{l[0]:>6} {l[1]:>5} {l[2]:>5}   1   {l[4]:>8}   0.0'
            new_top.append(newline)
            continue
            
        elif trigger == 'dihedrals':
            if int(l[4]) == 3:
                newline = f'{l[0]:>6} {l[1]:>5} {l[2]:>5} {l[3]:>5}   3     0.0   0.0   0.0   0.0   0.0   0.0'
                new_top.append(newline)
            elif int(l[4]) == 2:
                newline = f'{l[0]:>6} {l[1]:>5} {l[2]:>5} {l[3]:>5}   2   {l[5]:>8}     0.0'
                new_top.append(newline)
            elif (int(l[4]) == 1) or (int(l[4]) == 4):
                newline = f'{l[0]:>6} {l[1]:>5} {l[2]:>5} {l[3]:>5} {l[4]:>2}   {l[5]:>8}     0.0   {l[7]:>3}'
                new_top.append(newline)
            else:
                print(f'WARNING! The dihetral between atoms {l[0]} {l[1]} {l[2]} {l[3]} could not be set ot 0.0')
                print(f'The current force constant is {l[6]}; change its value to 0.0 in {name}.top, at line {i}')
            continue
        
        #if trigger == 'exclusions':
        #    new_top.append(line)

        elif trigger == 'pairs':
            newline = ";" + line
            new_top.append(newline)
        
        elif trigger == 'pairs_nb':
            newline = ";" + line
            new_top.append(newline)
  
        else:
            new_top.append(line)
                
    with open(f'topologies/{name}', "w") as f:
        for l in new_top:
            f.write(str(l) + '\n')


if __name__ == "__main__":
    feps, qmatoms, top, rs, ps, cutoff, prec = get_args()
    if prec not in ['single', 'double']:
        print("Precision can only be 'single' or 'double' (default: single).")
        sys.exit()

    qpdb, q1, q2, du, charges, soft_at, bevb, bconstr, soft_pairs, aevb, torevb, impevb = read_qm(qmatoms, rs, ps)

    if rs:
        try:
            rs_files = filelist(rs, top)
        except:
            print(f'{rs} files could not be opened')
            sys.exit()
    else:
        print("There are no parameter files for reactants state.")
        
    if ps:
        try:
            ps_files = filelist(ps)
        except:
            print(f'{ps} files could not be opened')
            sys.exit()
    else:
        print("There are no parameter files for products state.")

    try:
        atoms = atom_list(charges, qpdb)
    except:
        print('The van der Waals could not be read')
        sys.exit()

    try:
        # bonds_x2y contains the bonds that form or break
        bonds, bonds_x2y = bonds_list(rs, ps, rs_files, ps_files, qpdb, q1, q2, bevb)
    except:
        print('The bonds could not be read')
        sys.exit()

    try:
        # angles_x2y contains the angles that form or break
        angles, angles_x2y = angles_list(rs, ps, rs_files, ps_files, qpdb, q1, q2, aevb)
    except:
        print('No angles were found for this job')
        #sys.exit()

    try:
        torsions = torsions_list('torsions', rs, ps, rs_files, ps_files, qpdb, q1, q2, torevb)
    except:
        print('No torsions were found for this job')
    #    #sys.exit()

    try:
        impropers = torsions_list('impropers', rs, ps, rs_files, ps_files, qpdb, q1, q2, impevb)
    except:
        print('No impropers were found for this job')
        #sys.exit()

    ## pairs for 1-4 interaction and vdW params
    pairs_x2y, rs_vdw, ps_vdw = pairs_list(rs, ps, rs_files, ps_files, qpdb, q1, q2)

    # bonds and angles present in only one of the states
    bonds_solo, angles_solo, at_ad = coulomb_list(bonds_x2y, angles_x2y, pairs_x2y, bonds, angles)

    # soft core parameters
    try:
        soft, betas = soft_core(at_ad, bonds_solo, angles_solo, soft_at, soft_pairs)
    except:
        print('The soft-core parameters could not be generated')
        print('Check the input files and try again')
        sys.exit()

    # generate tabulated tables for soft-core potential
    gen_tables(betas, cutoff, prec)

    # excluded pairs due to the forming/breaking bonds
    nb55 = restore_nb(bonds_solo, angles_solo)

    # write topologies
    try:
        write_top(at_ad,top,atoms,bonds,angles,torsions,impropers,soft,bonds_solo,angles_solo,pairs_x2y,charges,rs_vdw,ps_vdw,feps,bconstr,nb55)
    except:
        print('Topology files could not be written')
        sys.exit()

    # write evbless.top
    try:
        evb_less('topol_000.top', 'evbless.top', du)
    except:
        print('evbless.top file could not be written')

    # show citing paper
    '''print("""
If you find these tools useful, please cite the following paper:
Gabriel Oanca, Florian van der Ent, Johan Åqvist, Efficient Empirical Valence Bond Simulations with GROMACS, Journal of Chemical Theory and Computation, 2023, doi: 10.1021/acs.jctc.3c00714
""")'''

