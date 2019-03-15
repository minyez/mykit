#!/usr/bin/env python3

import os
import subprocess as sp
import sys
from argparse import ArgumentParser

import numpy as np

from mykit.core.cell import sym_nat_from_atoms
from mykit.vasp.outcar import read_fermi
from mykit.vasp.poscar import poscar

# ===============================

def read_doscar(dosfile='DOSCAR'):

    # read in the DOSCAR file
    with open(dosfile,'r') as h_dos:
        doslines = h_dos.readlines()
    # split lines
    for i in range(len(doslines)):
        doslines[i] = doslines[i].split()

    # the first 5 lines include the number of atoms, the number of energy grids,
    # the coordinate system and the name of the system
    natoms_tot = int(doslines[0][0])
    nedos = int(doslines[5][2])
    # only co-linear spin is considered
    if len(doslines[6]) == 3:
        ispin = 1
    else:
        ispin = 2

    print("ISPIN = %d, EDOS = %d" % (ispin, nedos))

    energy_grid = [float(x[0]) for x in doslines[6:6+nedos]]
    energy_grid = np.array(energy_grid)

    # total DOS
    total_dos = []
    if ispin == 1:
    # total_dos is a single 1d numpy array for spin-unpolarized calculation
        total_dos = [float(x[1]) for x in doslines[6:6+nedos]]
    else:
    # total_dos is a nx2 numpy array for spin-polarized calculation
        for x in doslines[6:6+nedos]:
            total_dos.append([float(x[1]), float(x[2])]) 

    total_dos = np.array(total_dos)

    pdos = []

    for iatom in range(natoms_tot):
        pdos_atom = []
        for line in doslines[6+(iatom+1)*(nedos+1):6+(iatom+2)*(nedos+1)-1]:
            pdos_atom.append([float(x) for x in line[1:]])
        pdos_atom = np.array(pdos_atom)
        pdos.append(pdos_atom)

    return ispin, energy_grid, total_dos, pdos

# ===============================

def add_dos(atom_type_list, natoms_list, ispin, energy_grid, total_dos, pdos, \
            fermi_level=0.0, tot_only=False):

    l_dict = {0:'s',1:'p',2:'d',3:'f',4:'g',5:'h',6:'i'}

    ntype = len(natoms_list)
    nedos = len(energy_grid)
    atom_start_line = [sum(natoms_list[:iatom]) for iatom in range(ntype)]

    # set fermi level to energy zero
    energy_grid = energy_grid - fermi_level

    list_dos_sum_in_atom = []
    pdos_shape = pdos[0].shape

    # sum over atoms
    for itype in range(ntype):
        dos_sum_in_atom = np.zeros(pdos_shape, dtype='float64')
        for iatom in range(atom_start_line[itype], atom_start_line[itype]+natoms_list[itype]):
            dos_sum_in_atom = np.add(dos_sum_in_atom, pdos[iatom])
        list_dos_sum_in_atom.append(dos_sum_in_atom)

    if ispin == 1:
        lmax = int(np.sqrt(len(list_dos_sum_in_atom[0][0])) - 1)
    else:
        lmax = int(np.sqrt(len(list_dos_sum_in_atom[0][0])/2) - 1)
    print("lmax = ", lmax)

    # s-up, p-up, d-up, f-up, ..., s-dn, p-dn, d-dn, f-dn, ...
    list_dos_sum_in_atom_and_m = []

    # sum up according to the type of atom
    for itype in range(ntype):
        if ispin == 1:
            dos_sum_in_atom_and_m = []
            for l in range(lmax+1):
                temp_array = list_dos_sum_in_atom[itype][:,l*l:l*l+2*l+1].sum(axis=1)
                dos_sum_in_atom_and_m.append(temp_array)
            dos_sum_in_atom_and_m = np.array(dos_sum_in_atom_and_m).transpose()
            list_dos_sum_in_atom_and_m.append(dos_sum_in_atom_and_m)
        elif ispin == 2:
            dos_sum_in_atom_and_m = []
            for l in range(lmax+1):
                # spin up
                temp_array = list_dos_sum_in_atom[itype][:,2*l*l:2*l*l+(2*l+1)*2:2].sum(axis=1)
                dos_sum_in_atom_and_m.append(temp_array)
                # spin down
                temp_array = list_dos_sum_in_atom[itype][:,2*l*l+1:2*l*l+1+(2*l+1)*2:2].sum(axis=1)
                dos_sum_in_atom_and_m.append(temp_array)
            dos_sum_in_atom_and_m = np.array(dos_sum_in_atom_and_m).transpose()
            list_dos_sum_in_atom_and_m.append(dos_sum_in_atom_and_m)

    # export to dos.dat if ISPIN=1
    #        to dos.dat, dos_up.dat and dos_down if ISPIN=2, NOT IMPLEMENTED !

    if ispin == 1:
        with open('dos.dat','w') as h_dos:
            # write header
            h_dos.write("#energy    ")
            h_dos.write("        tot")
            if not tot_only:
                for itype in range(ntype):
                    for l in range(lmax+1):
                        h_dos.write("       %2s-%1s" % (atom_type_list[itype], l_dict[l]))
            h_dos.write("\n")

            for i in range(nedos):
                h_dos.write("%11.6f" % energy_grid[i])
                h_dos.write("%11.6f" % total_dos[i])
                if not tot_only:
                    for itype in range(ntype):
                        for l in range(lmax+1):
                            h_dos.write("%11.6f" % list_dos_sum_in_atom_and_m[itype][i][l])
                h_dos.write("\n")

# ===============================

def dos_anal_main(ArgList):

    description='''read DOSCAR to get total and atom-angular decomposed density of states (DOS).'''
    
    parser = ArgumentParser(description=description)
    parser.add_argument('--tot',dest='tot_only',action='store_true',help="flag to output total DOS only")
    # initialize options as 'opts'
    opts = parser.parse_args()
    

    fermi_level = read_fermi()

    atom_type_list, natoms_list = sym_nat_from_atoms(poscar.read_from_file().atoms)

    ispin, energy_grid, total_dos, pdos = read_doscar()

    add_dos(atom_type_list, natoms_list, ispin, energy_grid, total_dos, pdos, \
            fermi_level, tot_only=opts.tot_only)

# ===============================

if __name__ == "__main__":
    dos_anal_main(sys.argv)
