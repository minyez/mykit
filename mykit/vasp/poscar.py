# -*- coding: utf-8 -*-
'''Define the class for manipulating POSCAR, the VASP lattice input
'''
import logging
import os
import string

import numpy as np

from mykit.core.lattice import atoms_from_sym_nat, lattice, sym_nat_from_atoms
from mykit.core.utils import trim_comment


class PoscarError(Exception):
    pass

class poscar(lattice):
    '''The class to manipulate POSCAR, the VASP lattice input file.
    '''

    def __init__(self, cell, atoms, pos, **kwargs):
        super(poscar, self).__init__(cell, atoms, pos, **kwargs)

    # ? Rewrite atomType and typeIndex properties, as the element symbol can appear twice in POSCAR
    # ? But it may be good to let user deal with it, such as using "Fe1", "Fe2" to distinguish.
    # ? In this case, it should be careful to set POTCAR when recognizing atomic information in POSCAR

    def __print(self, fp):
            _cell, _atoms, _pos = self.get_latt()
            _syms, _nats = sym_nat_from_atoms(_atoms)
            print(self.comment, file=fp)
            print("1.00000", file=fp)
            for i in range(3):
                print("  %12.8f  %12.8f  %12.8f" \
                    % (self.cell[i,0], self.cell[i,1], self.cell[i,2]), file=fp)
            if _syms[0] != "A":
                print(*_syms, file=fp)
            print(*_nats, file=fp)
            if self.useSelDyn:
                print("Selective Dynamics", file=fp)
            print({"D": "Direct", "C": "Cart"}[self.coordSys], file=fp)
            for i in range(self.natoms):
                if self.useSelDyn:
                    _dyn = self.sdFlags(ia=i)
                else:
                    _dyn = []
                if _syms[0] != "A":
                    _ainfo = ['#{}'.format(_atoms[i])]
                else:
                    _ainfo = []
                _aflag = [{True:"T", False:"F"}[_d] for _d in _dyn] + _ainfo
                print("%15.9f  %15.9f %15.9f " % (self.pos[i,0], self.pos[i,1], self.pos[i,2]), \
                    *_aflag, file=fp)

    def print(self):
        '''Preview the POSCAR output
        '''
        from sys import stdout
        self.__print(stdout)

    def write(self, pathPoscar='POSCAR', backup=False, suffix="_bak"):
        '''Write POSCAR to path 
        '''
        # TODO test it
        _name = pathPoscar
        try:
            assert not os.path.isdir(_name)
        except AssertionError:
            raise PoscarError("The path to write POSCAR is a directory.")
        if os.path.isfile(_name) and backup:
            _bakname = _name + suffix.strip()
            os.rename(_name, _bakname)
        with open(_name, 'w') as f:
            self.__print(f)
    
    @classmethod
    def read_from_file(cls, pathPoscar="POSCAR"):
        '''Read from an existing POSCAR file ``pathPoscar``

        Args :
            pathPoscar (str) : path to the file to read as POSCAR
        '''
        # TODO may be decomposed to support reading structures stored in XDATCAR and OUTCAR?
        __flagSelDyn = False
        __fix = {}
        __fixDict = {'T': True, 'F': False}
        try:
            _f = open(pathPoscar, 'r')
        except FileNotFoundError as _err:
            raise PoscarError("Fail to open file: {}".format(pathPoscar))
        else:
            _symTypes = None
            # line 1: comment on system
            _comment = _f.readline().strip()
            # line 2: scale
            _scale = float(_f.readline().strip())
            # line 3-5: lattice vector
            _cell = [_f.readline().split() for i in range(3)]
            try:
                _cell = np.array(_cell, dtype=cls._dtype) * _scale
            except ValueError:
                _f.close()
                raise PoscarError("Bad lattice vector: {}".format(pathPoscar))
            # Next 2 or 1 line(s), depend on whether element symbols are typed or not
            _line = _f.readline().strip()
            if _line[0] in string.ascii_letters:
                _symTypes = _line.split()
                _line = _f.readline().strip()
            if _line[0] in string.digits[1:]:
                _natomsType = [int(_x) for _x in _line.split()]
                if _symTypes is None:
                    cls.print_cm_warn("No atom information in POSCAR: {}".format(pathPoscar))
                    _symTypes = [string.ascii_lowercase[_i] for _i,_x in enumerate(_natomsType)]
            else:
                _f.close()
                raise PoscarError("Bad POSCAR atomic format: {}".format(pathPoscar))
            try:
                assert len(_symTypes) == len(_natomsType)
            except AssertionError:
                _f.close()
                raise PoscarError("Inconsistent input of symbol and numbers of atom: {}".format(pathPoscar))
            _atoms = atoms_from_sym_nat(_symTypes, _natomsType)
            _natoms = sum(_natomsType)
            # Next 2 or 1 line(s), depend on whether 'selective dynamics line' is typed
            _line = _f.readline().strip()
            if _line[0].upper() == "S":
                __flagSelDyn = True
                _line = _f.readline().strip()
            if _line[0].upper() in ["C", "K", "D"]:
                _cs = {"C":"C", "K":"C", "D":"D"}[_line[0].upper()]
            else:
                raise PoscarError("Bad coordinate system: {}".format(pathPoscar))
            # Next _natoms lines: read atomic position and selective dynamics flag
            _pos = []
            _mult = 1.0E0
            if _cs == "C":
                _mult = _scale
            for _i in range(_natoms):
                try:
                    _line = _f.readline().strip()
                    _words = trim_comment(_line, r'[\#]').split()
                    _pos.append([float(_x) for _x in _words[:3]])
                    if __flagSelDyn:
                        __fixFlag = [__fixDict.get(_words[i]) for i in range(3,6)]
                        if None in __fixFlag:
                            raise IndexError
                        elif __fixFlag == [True, True, True]:
                            pass
                        else:
                            __fix.update({_i: __fixFlag})
                except ValueError:
                    _f.close()
                    raise PoscarError("Bad internal coordinates at atom line {}: {}".format(_i+1, pathPoscar))
                except IndexError:
                    _f.close()
                    raise PoscarError("Bad selective dynamics flag at atom line {}: {}".format(_i+1, pathPoscar))
            _pos = np.array(_pos, dtype=cls._dtype) * _mult
            _f.close()
            return cls(_cell, _atoms, _pos, unit="ang", coordSys=_cs, allRelax=True, selectDyn=__fix, comment=_comment)

    @classmethod
    def create_from_lattice(cls, latt):
        '''Create POSCAR from ``lattice`` instance ``latt``.
        '''
        try:
            assert isinstance(latt, lattice)
        except AssertionError:
            raise PoscarError("the input is not a lattice instance")
        __kw = latt.get_kwargs()
        return cls(*latt.get_latt(), **__kw)

    # TODO inherit factory methods from ``lattice`` for space group and Bravis cell creation
