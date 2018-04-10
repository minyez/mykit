#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pc_readce.py
# Creation Date : 09-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function
import sys
import string

Periodic_Table = (
                'H' , 'He', 
                'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', 
                'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 
                'K' , 'Ca', 
                'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 
                'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
                'Rb', 'Sr', 
                'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
                'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe', 
                'Cs', 'Ba', 
                'La', 
                'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 
                'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 
                'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 
                'Fr', 'Ra', 
                'Ac', 
                'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 
                'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 
                'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' 
              )


def common_read_chemical_formula(chem_formula, debug=False):
    
    atom_type = []
    natom_list = []


    length_symbol_tmp = 0
    natom = 0
    symbol_tmp = ''
    str_natom_tmp = ''

    # procedure: 
    for i in range(len(chem_formula)):
        char_ce = chem_formula[i]
        if char_ce in string.letters:
            # If an uppercase is met, save the last element if str_ce is symbol_tmp is not empty
            # and start a new element
            if char_ce in string.uppercase:
                if len(symbol_tmp) > 0:
                    atom_type.append(symbol_tmp)
                    natom_list.append(int(str_natom_tmp))
                    # clear the natom and length counter
                    length_symbol_tmp = 0
                    str_natom_tmp = ''

                symbol_tmp = char_ce

            elif char_ce in string.lowercase:
                symbol_tmp += char_ce
            length_symbol_tmp += 1
            # no chemical symbol has a length larger than 2
            assert length_symbol_tmp <= 2
        elif char_ce in string.digits:
            str_natom_tmp += char_ce
            
        # save the last element
        if i == (len(chem_formula)-1):
            atom_type.append(symbol_tmp)
            natom_list.append(int(str_natom_tmp))

    if debug:
        print(atom_type, natom_list)

    return atom_type, natom_list
            

#
## ==============================
#
#if __name__ == "__main__":
#    if len(sys.argv) > 3:
#        print("Error: too many argumetns. Need the chemical formula only, e.g. Si1O2")
#        sys.exit(1)
#    elif len(sys.argv) == 3:
#        common_read_chemical_formula(sys.argv[1], debug=sys.argv[2])
#    elif len(sys.argv) == 2:
#        common_read_chemical_formula(sys.argv[1])
#    elif len(sys.argv) == 1:
#        # print a simple help
#        print("Usage: pc_readce.py chem_formula [debug]")
#        sys.exit(2)
#
