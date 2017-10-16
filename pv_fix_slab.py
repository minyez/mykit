#!/usr/bin/env python
#
#

# * * * * * * * * * * * * * * * * * * *
#
#     File Name :
# Creation Date :
# Last Modified : Mon 16 Oct 2017 03:54:59 PM CST
#       Purpose : sort atoms and fix some layers according to needs,
#                 modified from the script written by Yue-Chao Wang
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# * * * * * * * * * * * * * * * * * * *

import sys
import argparse

# ====================================================

def get_atom_info(poscar_lines):

    Atom_Type = poscar_lines[5].split()
    Atom_Numb = [int(x) for x in poscar_lines[6].split()]
    Atom_St_Line = [8]
    for i in xrange(len(Atom_Type)-1):
        Atom_St_Line.append(Atom_St_Line[i]+Atom_Numb[i])
#    print Atom_Type
#    print Atom_Numb
#    print Atom_St_Line

    return  Atom_Type, Atom_Numb, Atom_St_Line

# ====================================================

def sort_xyz(poscar_lines,direct,debug=False):
    '''
    Sort atoms of each type
    '''
    if poscar_lines[7].split()[0].startswith('S'):
        del poscar_lines[7]
    Atom_Type, Atom_Numb, Atom_St_Line = get_atom_info(poscar_lines)

    for i in xrange(len(Atom_Type)):
        Coords = [float(x.split()[direct]) for x in poscar_lines[Atom_St_Line[i]:Atom_St_Line[i]+Atom_Numb[i]]]
        Index = sorted(range(len(Coords)),key=lambda k:Coords[k])
        if debug: print Index
        temp_list = []

# sort the coordinates for each type of atoms
        for j in xrange(len(Coords)):
            temp_list.append(poscar_lines[Atom_St_Line[i]+Index[j]])
        if debug: print temp_list

        poscar_lines[Atom_St_Line[i]:Atom_St_Line[i]+Atom_Numb[i]] = temp_list

    Coords_all = [float(x.split()[direct]) for x in poscar_lines[Atom_St_Line[0]:Atom_St_Line[-1]+Atom_Numb[-1]]]
    all_atom_sort = sorted(range(len(Coords_all)),key=lambda k:Coords_all[k])
    if debug:
        print Atom_St_Line
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

def fix_xyz_range(poscar_lines,start,end,direct,reverse=False,debug=False):
    '''
    Fix x/y/z between start and end points. Need sorted poscar_lines
    Use reverse if you want to the points out of this range to be fixed. which is rarely used
    '''
    poscar_lines[7] = 'Selective dynamics\n'+poscar_lines[7]

    TrSurf, FlsSurf = flag_fix(reverse)
    Atom_Type, Atom_Numb, Atom_St_Line = get_atom_info(poscar_lines)

    n_line = 8
    for i in xrange(len(Atom_Numb)):
        while n_line < Atom_St_Line[i]+Atom_Numb[i]:
            Coord_line = poscar_lines[n_line].split()
            key = float(Coord_line[direct])
            # in the range
            if key >= start and key <= end:
                poscar_lines[n_line] = '   '+'   '.join(poscar_lines[n_line].split()[0:3])+FlsSurf+'# %s\n' % Atom_Type[i]
            # out of the range
            else:
                poscar_lines[n_line] = '   '+'   '.join(poscar_lines[n_line].split()[0:3])+TrSurf+'# %s\n' % Atom_Type[i]
            n_line = n_line + 1

# ====================================================

def fix_xyz_number_bottom(poscar_lines,sorted_list,fixnum,reverse=False,debug=False):
    '''
    Fix atoms from bottom of slab. Need sorted poscar_lines.
    Sorted list
    Use reverse if you want to the points out of this range to be fixed. which is rarely used
    '''
    poscar_lines[7] = 'Selective dynamics\n'+poscar_lines[7]

    TrSurf, FlsSurf = flag_fix(reverse)
    Atom_Type, Atom_Numb, Atom_St_Line = get_atom_info(poscar_lines)
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
        while n_line < Atom_St_Line[i]+Atom_Numb[i]:
            if (n_line - 8) in fix_index:
                poscar_lines[n_line] = '   '+'  '.join(poscar_lines[n_line].split()[0:3])+FlsSurf+'# %s\n' % Atom_Type[i]
            else:
                poscar_lines[n_line] = '   '+'  '.join(poscar_lines[n_line].split()[0:3])+TrSurf+'# %s\n' % Atom_Type[i]
            n_line = n_line + 1

# ====================================================

def fix_xyx_number_sym(poscar_lines,sorted_list,surf_num,reverse=True,debug=False):
    '''
    Fix atoms symmetrically.
    surf_num is the number of atoms to relax in each of surface
    '''
    poscar_lines[7] = 'Selective dynamics\n'+poscar_lines[7]

    TrSurf, FlsSurf = flag_fix(reverse)
    Atom_Type, Atom_Numb, Atom_St_Line = get_atom_info(poscar_lines)

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
        while n_line < Atom_St_Line[i]+Atom_Numb[i]:
            if (n_line - 8) in fix_index:
                poscar_lines[n_line] = '   '+'  '.join(poscar_lines[n_line].split()[0:3])+FlsSurf+'# %s\n' % Atom_Type[i]
            else:
                poscar_lines[n_line] = '   '+'  '.join(poscar_lines[n_line].split()[0:3])+TrSurf+'# %s\n' % Atom_Type[i]
            n_line = n_line + 1

# ====================================================

def Main(ArgList):

# =================== Parser ==========================

    parser = argparse.ArgumentParser(description="Fix Coordinates")
    group1 = parser.add_mutually_exclusive_group()

    parser.add_argument("-f",dest="filename",default="POSCAR",help="File Name Default=POSCAR")
    group1.add_argument("-n",dest="fixnum",type=int,help="Number of atoms from bottom to set F/T")
    group1.add_argument("-r",dest="range",type=float,nargs=2,help="The range of F or T")
    parser.add_argument("-s",dest="sym",help="Flag for symmetric fixing.",action='store_true')
    parser.add_argument("-d",dest="direction",default="z",help="Sort according to which direction x, y or z")
    parser.add_argument("-t",action='store_true',help="Using this Tag Set T in Range, Default is False")
    parser.add_argument("-o",dest="Output",default="POSCAR_new",help="Name of New POSCAR")
#    parser.add_argument("-l",action="store_true",help="List Coordinates after sorted")
    parser.add_argument("-D",action="store_true",help="Debug mode")

    opts = parser.parse_args()

# ====================================================

    NewFileName = opts.Output
    FileName = opts.filename

    PosFile = open(FileName,"r")
    ofile = open(NewFileName,"w+")
    lines = PosFile.readlines()
    PosFile.close()

    if opts.direction == 'x':
        drct = 0
    elif opts.direction == 'y':
        drct = 1
    elif opts.direction == 'z':
        drct = 2
    else :
        print "Wrong in PUT"
        exit()

    all_atom_sort = sort_xyz(lines,drct,opts.D)
    if opts.range is not None:
        start = min(opts.range)
        end = max(opts.range)
        fix_xyz_range(lines,start,end,drct,opts.t,opts.D)
    elif opts.fixnum is not None:
        if not opts.sym:
            fix_xyz_number_bottom(lines,all_atom_sort,opts.fixnum,opts.t,opts.D)
        else:
            print "Symmetric fixing ..."
            fix_xyx_number_sym(lines,all_atom_sort,opts.fixnum,True,opts.D)
    else:
        print "Specify either range of number of atoms to fix"
        exit()

    for i in xrange(len(lines)):
        ofile.write(lines[i])

    ofile.close()

# ====================================================

if __name__ == "__main__":
    Main(sys.argv)
