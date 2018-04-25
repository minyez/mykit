#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pc_elements.py
# Creation Date : 09-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function
import sys
import string

# ====================================================

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

# ====================================================
# Standard atomic weight, or relative atomic mass of the element 
# Data from NIST database 144, last update January 2015
# https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
# The data is double-checked with Wolfram Language ElementData (WLE)
# In case that an atomic-weight interval is met in NIST, the data in WLE is used instead and noted below.
# All data are accurate to 6 decimals, if available
# WLE data used: H, Li, B, C, N, O, Mg, Si, S, Cl, Br, Tl
Atom_Weight = (
1.008     , 4.002602 ,
6.94      , 9.012183 , 10.81    , 12.011 , 14.007   , 15.999 , 18.998403, 20.1797   ,
22.989769 , 24.305   , 26.981539, 28.085 , 30.973762, 32.06  , 35.45    , 39.948    ,
39.0983   , 40.078   ,
44.955908 , 47.867   , 50.9415  , 51.9961, 54.938044, 55.845 , 58.933194, 58.6934   , 63.546  , 65.38    ,
69.723    , 72.6308  , 74.921595, 78.971 , 79.904   , 83.798 ,
85.4678   , 87.62    ,
88.90584  , 91.224   , 92.90637 , 95.95  , 98       , 101.07 , 102.90550, 106.42    , 107.8682, 112.414  ,
114.818   , 118.710  , 121.760  , 127.60 , 126.90447, 131.293,
132.905452, 137.327  ,
138.90547 ,
140.116   , 140.90766, 144.242  , 145    , 150.36   , 151.964, 157.25   , 158.92535 , 162.500 , 164.93033, 167.259, 168.93422, 173.054, 174.9668,
178.49    , 180.94788, 183.84   , 186.207, 190.23   , 192.217, 195.084  , 196.966569, 200.592 ,
204.38    , 207.2    , 208.98040, 209    , 210      , 222    ,
223       , 226      ,
227       ,
232.0377  , 231.03588, 238.02891, 237    , 244 
)   # End at Pu

# ====================================================

def common_molar_mass_calculator(chem_formula):

    atom_type, natom_list = common_read_chemical_formula(chem_formula, debug=False)

    molar_mass = 0.0

    for i in range(len(natom_list)):
        atom_num = Periodic_Table.index(atom_type[i])
        molar_mass += Atom_Weight[atom_num] * natom_list[i]
   
    return molar_mass

# ====================================================

def common_read_chemical_formula(chem_formula, debug=False):
    
    atom_type = []
    natom_list = []

    length_symbol_tmp = 0
    natom = 0
    symbol_tmp = ''
    str_natom_tmp = ''

    for i in range(len(chem_formula)):
        char_ce = chem_formula[i]
        if char_ce in string.letters:
            # If an uppercase is met, save the last element if symbol_tmp is not empty
            # and start a new element
            if char_ce in string.uppercase:
                if len(symbol_tmp) > 0:
                    # check if the symbol is a valid element symbol
                    if symbol_tmp in Periodic_Table:
                        atom_type.append(symbol_tmp)
                        natom_list.append(int(str_natom_tmp))
                    else:
                        raise ValueError("Invalid element symbol not found in the periodic table.")
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
        else:
            raise TypeError("Invalid character for chemical formula.")
            
        # save the last element
        if i == (len(chem_formula)-1):
            atom_type.append(symbol_tmp)
            natom_list.append(int(str_natom_tmp))

    if debug:
        print(atom_type, natom_list)

    return atom_type, natom_list
