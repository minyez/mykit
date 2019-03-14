#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pc_crystal.py
# Creation Date : 20-06-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

from __future__ import print_function
import sys,os
from pv_classes import vasp_read_poscar
from argparse import ArgumentParser
import numpy as np
try:
    import spglib
except ImportError:
    import pyspglib as spglib

def __check_input_type(filename):

    intype = None

    for t_vasp in ['.vasp', 'poscar']:
        if filename.lower().endswith(t_vasp):
            intype = 'vasp'
            break

    for t_w2k in ['.struct']:
        if filename.lower().endswith(t_w2k):
            intype = 'w2k'
            break

    return intype


def __cell_from_vasp(filename):

    poscar = vasp_read_poscar(filename, verbose=False)
    poscar.action_cart2dirt()
    lattice = poscar.lattice
    pos = np.array(poscar.innerpos)

    nums = []
    for i in range(len(poscar.atom_num)):
        nums_i = [i,] * poscar.atom_num[i]
        nums.extend(nums_i)

    return (lattice, pos, nums)


def __cell_from_w2k(filename):
    raise IOError("WIEN2k input has not yet been supported. Exit. ")

def __cell_from_intype(intype, filename):

    cell = None

    if intype.lower() in ['v', 'vasp']:
        cell = __cell_from_vasp(filename)
    if intype.lower() in ['w', 'w2k', 'wien2k', 'wien']:
        cell = __cell_from_w2k(filename)

    return cell 

def __print_sym_info(cell, print_std, print_symop):
    '''
    Return symmtery information of the cell tuple
    '''
    print("===== Crystall Symmetry =====")
    ds = spglib.get_symmetry_dataset(cell)
    print("Space group: %8s" % ds['international'])
    print("     Number: %8d" % ds['number'])
    #print("Point group: %10s" % spglib.get_pointgroup(cell))
    print("Point group: %8s" % ds['pointgroup'])
    origin = ds['origin_shift']
    print("     Origin: [%5.3f, %5.3f, %5.3f]" % tuple(origin))
    rot = ds['rotations']
    tau = ds['translations']
    # symmtry operations
    print("# Sym. Ops.: %8d" % len(rot))
    if print_symop:
        sym_in_line = 6
        print("  Sym. Ops.: (%d in a line)" % sym_in_line)
        
        left = len(rot) % sym_in_line
        line = (len(rot)-left) / sym_in_line
        for i in range(line):
            print("="*(9+6+1)*sym_in_line)
            for j in range(3):
                for k in range(sym_in_line):
                    print("%3d%3d%3d" % (rot[i*sym_in_line+k][j][0], rot[i*sym_in_line+k][j][1], rot[i*sym_in_line+k][j][2]) ,end='')
                    print("|%4.2f " % tau[i*sym_in_line+k][j] ,end='')
                    if j == 1:
                        print(",", end='')
                    else:
                        print(" ", end='')
                print("")
        print("="*(9+6+1)*sym_in_line)

        if left != 0:
            for j in range(3):
                for i in range(left):
                    sym_index = i - left
                    print("%3d%3d%3d" % (rot[sym_index][j][0], rot[sym_index][j][1], rot[sym_index][j][2]) ,end='')
                    print("|%4.2f " % tau[sym_index][j], end='')
                    if j == 1 and i != left - 1:
                        print(",", end='')
                    else:
                        print(" ", end='')
                print("")

    if print_std:
        print(ds['std_lattice'])
        print(ds['std_positions'])
        print(ds['std_types'])
    if not print_symop:
        print("=============================")


def __Main(ArgList):

    description='''Return symmtery information of crystal structure by taking advantage of spglib. Currently support VASP only.'''
    
    parser = ArgumentParser(description=description)
    #group1 = parser.add_mutually_exclusive_group()
    
    parser.add_argument("filename", nargs=1, help="input file containing crystal structure")
    parser.add_argument("-i", dest="intype", default=None, help="type of crystal input. default automatic detect")
    parser.add_argument("-D", dest="debug", action="store_true", help="debug mode")
    parser.add_argument("--std", dest="std", action="store_true", help="flag to return standardized cell info")
    parser.add_argument("--sym", dest="sym", action="store_true", help="flag to show the symmtry operations")
    
    # initialize options as 'opts'
    opts = parser.parse_args()
    filename = opts.filename[0]

    if not os.path.exists(filename):
        raise IOError('Input file %s does not exist.' % filename)

    if opts.intype is None:
        intype = __check_input_type(filename)
    else:
        intype = opts.intype

    cell =  __cell_from_intype(intype, filename)
    if opts.debug:
        print(cell)
    # space group
    __print_sym_info(cell, opts.std, opts.sym)


# ==============================

if __name__ == "__main__":
    __Main(sys.argv)

