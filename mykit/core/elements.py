# coding=utf-8

from __future__ import print_function
import string

# pylint: disable=bad-continuation,bad-whitespace,line-too-long
periodicTable = (
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
# In case that an atomic-weight interval is met in NIST, 
# the data in WLE is used instead and noted below.
# All data are accurate to 6 decimals, if available
# WLE data used: H, Li, B, C, N, O, Mg, Si, S, Cl, Br, Tl
atomWeight = (
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

def common_molar_mass_calculator(chemFormula):
    '''calculate molar weight for a chemical formula

    Args:
        chemFormula (str) : Chemical formula of compound. Parentheses are not supported.
    '''

    _at, _natomList, _compo = common_read_chemical_formula(chemFormula, debug=False)

    _molarMass = 0.0

    for _i, _nat in enumerate(_natomList):
        _atomNum = periodicTable.index(_at[_i])
        _molarMass += atomWeight[_atomNum] * _nat

    return _molarMass

# ====================================================

def common_read_chemical_formula(chemFormula, debug=False):

    '''Read the chemical formula

    Args:
        chemFormula (str) : the chemical formula of the system which reflects the total numbers of each element
            in the molecule or the unit cell of crystal.
        debug (bool) : the flag to switch on debug mode.
        
    Returns:
        atomType (list of str) : the types of elements in the system.
        natomList (list of int) : the numbers of atoms for each type of element correspondent to the atomType.
        compositions (int) : the number of compositions in the system.
    '''
    
    atomType = []
    natomList = []

    _lenSymbol = 0
    _natom = 0
    _symbol = ''
    _strNatom = ''

    for i, _charEle in enumerate(chemFormula):
        if _charEle in string.ascii_letters:
            # If an uppercase is met, save the last element if _symbol is not empty
            # and start a new element
            if _charEle in string.ascii_uppercase:
                if len(_symbol) > 0:
                    # check if the symbol is a valid element symbol
                    if _symbol in periodicTable:
                        atomType.append(_symbol)
                        natomList.append(int(_strNatom))
                    else:
                        raise ValueError("Invalid element symbol not found in the periodic table.")
                    # clear the natom and length counter
                    _lenSymbol = 0
                    _strNatom = ''

                _symbol = _charEle

            elif _charEle in string.ascii_lowercase:
                _symbol += _charEle
            _lenSymbol += 1
            # no chemical symbol has a length larger than 2
            assert _lenSymbol <= 2
        elif _charEle in string.digits:
            _strNatom += _charEle
        else:
            raise TypeError("Invalid character for chemical formula.")
            
        # save the last element
        if i == (len(chemFormula)-1):
            atomType.append(_symbol)
            natomList.append(int(_strNatom))

    compositions = len(atomType) 
        
    if debug:
        print(atomType, natomList)

    return atomType, natomList, compositions
