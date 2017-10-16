#!/usr/bin/env python
#
# Read XDATCAR and get intermediate structures during optimization
#
# 2017.05.09

import os
from pv_calc_utils import common_io_cleandir

struct_dir = 'inter_structs'

common_io_cleandir(struct_dir)
os.chdir(struct_dir)
with open('../XDATCAR','r') as f:
    lines = f.readlines()
    Atoms_Num = lines[6].split()
    head = ''.join(lines[0:7])+'Direct\n'
    natoms = sum(int(x) for x in Atoms_Num)
#    print head
#    print natoms
    n_config = 0
    nline = 8
    while nline < len(lines):
        with open('POSCAR_'+str(n_config)+'.vasp','w') as ofile:
            ofile.write(head)
            Coord = ''.join(lines[nline:nline+natoms])
            ofile.write(Coord)
        nline = nline + natoms + 1
        n_config = n_config + 1
#        print nline,n_config
