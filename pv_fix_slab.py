#!/usr/bin/env python
#
# ====================================================
#
#     File Name : pv_fix_slab.py
# Creation Date : 2017-06-20
# Last Modified : Mon 30 Oct 2017 09:31:28 PM CST
#       Purpose : sort atoms and fix some layers according to needs,
#                 modified from the script written by Yue-Chao Wang
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

import sys
import argparse
from pv_classes import vasp_read_poscar
from pv_slab_utils import sort_coor, flag_fix

# ====================================================

def get_atom_info(poscar):
    '''
    read in the atomic information from POSCAR
    there should be 8 lines in the heading:
    1 for comment, 1 for scale, 3 for lattice vectors
    2 for atoms type and numbers
    '''
    Atom_Type = poscar.atom_type
    Atom_Numb = poscar.atom_num
    Atom_Line = [8]
    for i in xrange(len(Atom_Type)-1):
        Atom_Line.append(Atom_Line[i]+Atom_Numb[i])
    return  Atom_Type, Atom_Numb, Atom_Line

# ====================================================

def fix_xyz_range(poscar,start,end,direct,reverse=False,debug=False):
    '''
    Fix x/y/z between start and end points. Need sorted poscar.poslines
    Use reverse if you want to the points out of this range to be fixed. which is rarely used
    '''

    idir = direct - 1
    TrSurf, FlsSurf = flag_fix(reverse)
    Atom_Type, Atom_Numb, Atom_Line = get_atom_info(poscar)

    rewrite_poscar(poscar,Atom_Numb,Atom_Type,fix_index,TrSurf,FlsSurf)
    n_line = 8
    for i in xrange(len(Atom_Numb)):
        while n_line < Atom_Line[i]+Atom_Numb[i]:
            Coord_line = poscar.poslines[n_line].split()
            key = float(Coord_line[idir])
            # in the range
            if key >= start and key <= end:
                rewrite_posline(poscar,n_line,FlsSurf,Atom_Type[i])
            # out of the range
            else:
                rewrite_posline(poscar,n_line,TrSurf,Atom_Type[i])
            n_line = n_line + 1

# ====================================================

def fix_xyz_number_bottom(poscar,sorted_list,fixnum,reverse=False,debug=False):
    '''
    Fix atoms from bottom of slab. Need sorted poscar.poslines.
    Sorted list
    Use reverse if you want to the points out of this range to be fixed. which is rarely used
    '''

    TrSurf, FlsSurf = flag_fix(reverse)
    Atom_Type, Atom_Numb, Atom_Line = get_atom_info(poscar)
    if debug:
        print sorted_list
    if fixnum <= 0:
        fix_index = []
    elif fixnum <= len(sorted_list):
        fix_index = sorted_list[0:fixnum]
    else:
        fix_index = sorted_list

    if reverse:
        print "Fixing from the bottom..."
    else:
        print "Fixing from the top..."

    n_line = 8
    for i in xrange(len(Atom_Numb)):
        while n_line < Atom_Line[i]+Atom_Numb[i]:
            if (n_line - 8) in fix_index:
                rewrite_posline(poscar,n_line,FlsSurf,Atom_Type[i])
            else:
                rewrite_posline(poscar,n_line,TrSurf,Atom_Type[i])
            n_line = n_line + 1


# ====================================================

def fix_xyx_number_sym(poscar,sorted_list,surf_num,reverse=True,debug=False):
    '''
    Fix atoms symmetrically.
    surf_num is the number of atoms to relax in each of surface
    '''
    TrSurf, FlsSurf = flag_fix(reverse)
    Atom_Type, Atom_Numb, Atom_Line = get_atom_info(poscar)

    if debug:
        print sorted_list
    if surf_num <= 0:
        fix_index = []
    elif 2*surf_num <= len(sorted_list):
        fix_index = sorted_list[0:surf_num] + sorted_list[-surf_num:]
    else:
        fix_index = sorted_list

    n_line = 8
    for i in xrange(len(Atom_Numb)):
        while n_line < Atom_Line[i]+Atom_Numb[i]:
            if (n_line - 8) in fix_index:
                rewrite_posline(poscar,n_line,FlsSurf,Atom_Type[i])
            else:
                rewrite_posline(poscar,n_line,TrSurf,Atom_Type[i])
            n_line = n_line + 1

def rewrite_posline(poscar,n_line,FlagSurf,atomtype):
    poscar.poslines[n_line] =  "%17.8f%17.8f%17.8f" % \
                    (poscar.innerpos[n_line-8][0],poscar.innerpos[n_line-8][1],poscar.innerpos[n_line-8][2]) \
                    + FlagSurf + '#%2s\n' % atomtype

# ====================================================

def Main(ArgList):

# =================== Parser ==========================

    parser = argparse.ArgumentParser(description="Fix Coordinates")
    group1 = parser.add_mutually_exclusive_group()

    parser.add_argument("-f",dest="filename",default="POSCAR",help="File Name Default=POSCAR")
    group1.add_argument("-n",dest="fixnum",type=int,help="Number of atoms from bottom to set F/T")
    group1.add_argument("-r",dest="range",type=float,nargs=2,help="The range of F or T")
    parser.add_argument("-s",dest="sym",help="Flag for symmetric fixing.",action='store_true')
    parser.add_argument("-d",dest="direction",type=int,default="3",help="Sort according to which direction 1, 2 or 3")
    parser.add_argument("-t",action='store_true',help="Using this Tag Set T in Range, Default is False")
    parser.add_argument("-o",dest="Output",default="POSCAR_new",help="Name of New POSCAR")
    parser.add_argument("-D",action="store_true",help="Debug mode")

    opts = parser.parse_args()

# ====================================================

    NewFileName = opts.Output
    FileName = opts.filename

    poscar = vasp_read_poscar(FileName)

    drct = opts.direction
    if drct not in [1,2,3]:
        print "Wrong in PUT"
        sys.exit(1)

# all_atom_sort is a sorted list
# the poscar.poslines is already sorted after the sort_coor function
#    all_atom_sort = sort_coor(poscar,drct)
    all_atom_sort = poscar.action_sort_coor(drct)
    if opts.range is not None:
        start = min(opts.range)
        end = max(opts.range)
        fix_xyz_range(poscar,start,end,drct,opts.t,opts.D)
    elif opts.fixnum is not None:
        if not opts.sym:
            fix_xyz_number_bottom(poscar,all_atom_sort,opts.fixnum,opts.t,opts.D)
        else:
            print "Symmetric fixing ..."
            fix_xyx_number_sym(poscar,all_atom_sort,opts.fixnum,True,opts.D)
    else:
        print "Specify either range of number of atoms to fix"
        exit()

    poscar.write_poscar(opts.Output)

# ====================================================

if __name__ == "__main__":
    Main(sys.argv)
