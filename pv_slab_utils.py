#!/usr/bin/env python
#
# ====================================================

import sys, os, shutil, subprocess, copy
import commands, string
import numpy as np
from argparse import ArgumentParser

# ====================================================

def vasp_get_NPAR(np):

    npar = 1
    nsqrt_max = int(sqrt(np)) + 1
    nsqrt = nsqrt_max
    while nsqrt > 0:
        if nsqrt % 2 == 0 and np % nsqrt == 0:
            npar = nsqrt
            break
        nsqrt = nsqrt - 1
    return npar

# ====================================================

def sort_coor(poscar,direct,debug=False):
    '''
    Sort atoms of each type
    poscar is the 'vasp_read_poscar' class
    '''
    Atom_Type = poscar.atom_type
    Atom_Numb = poscar.atom_num
    # skip the heading
    Atom_Line = [8]
    for i in xrange(len(Atom_Type)-1):
        Atom_Line.append(Atom_Line[i]+Atom_Numb[i])

    for i in xrange(len(Atom_Type)):
        Coords = [float(x.split()[direct]) for x in poscar.poslines[Atom_Line[i]:Atom_Line[i]+Atom_Numb[i]]]
        Index = sorted(range(len(Coords)),key=lambda k:Coords[k])
        if debug: print Index
        temp_list = []

# sort the coordinates for each type of atoms
        for j in xrange(len(Coords)):
            temp_list.append(poscar.poslines[Atom_Line[i]+Index[j]])
        if debug: print temp_list

        poscar.poslines[Atom_Line[i]:Atom_Line[i]+Atom_Numb[i]] = temp_list

    poscar.write_innerpos_from_lines()
    Coords_all = [float(x.split()[direct]) for x in poscar.poslines[Atom_Line[0]:Atom_Line[-1]+Atom_Numb[-1]]]
    all_atom_sort = sorted(range(len(Coords_all)),key=lambda k:Coords_all[k])
    if debug:
        print Atom_Line
        print len(Coords_all)
        print all_atom_sort
        print Coords_all

    return all_atom_sort

# ====================================================

def flag_fix(reverse=False):

    FlsSurf = "  F  F  F "
    TrSurf = "  T  T  T "
    if reverse:
        TrSurf = "  F  F  F "
        FlsSurf = "  T  T  T "
    return TrSurf, FlsSurf

# ====================================================

