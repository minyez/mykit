# -*- coding: utf-8 -*-
'''Define the class for manipulating POSCAR, the VASP lattice input
'''
import numpy as np
import logging
import string
from mykit.core.lattice import lattice, atoms_from_sym_nat

class poscarError(Exception):
    pass

class poscar(lattice):
    '''The class to manipulate POSCAR, the VASP lattice input file.
    '''

    def __init__(self, cell, atoms, pos, **kwargs):
        super(poscar, self).__init__(cell, atoms, pos, **kwargs)

    #TODO rewrite atomType and typeIndex properties
    
    @classmethod
    def read_from_file(cls, poscarPath):
        '''Read from an existing POSCAR file ``poscarPath``

        Args :
            poscarPath (str) : path to the file to read as POSCAR
        '''
        # TODO can be decomposed to support reading structures stored in XDATCAR and OUTCAR
        try:
            _f = open(poscarPath, 'r')
        except FileNotFoundError as _err:
            raise poscarError("Fail to open file: {}".format(poscarPath))
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
                raise poscarError("Bad lattice vector")
            # Next 2 or 1 line(s), depend on whether element symbols are typed or not
            _line = _f.readline().strip()
            if _line[0] in string.ascii_letters:
                _symTypes = _line.split()
                _line = _f.readline()
            if _line[0] in string.digits[1:]:
                _natomsType = [int(_x) for _x in _line.split()]
                if _symTypes is None:
                    cls.print_cm_warn("No atom information in POSCAR")
                    _symTypes = [string.ascii_lowercase[_i] for _i,_x in enumerate(_natomsType)]
            else:
                _f.close()
                raise poscarError("Bad POSCAR format: {}".format(poscarPath))
            try:
                assert len(_symTypes) == len(_natomsType)
            except AssertionError:
                _f.close()
                raise poscarError("Inconsistent input of symbol and numbers of atom")
            _atoms = atoms_from_sym_nat(_symTypes, _natomsType)
            _natoms = sum(_natomsType)
            # Next 2 or 1 line(s), depend on whether 'selective dynamics line' is typed
            _fSelectDyn = False
            _line = _f.readline().strip()
            if _line[0].upper() not in ["C", "D", "K"]:
                _fSelectDyn = True
                _line = _f.readline().strip()
            if _line[0].upper() in ["C", "K"]:
                _cs = "C"
            else:
                _cs = "D"
            # Next _natoms lines: read atomic position
            # TODO extract selective dynamics information
            _pos = []
            _mult = 1.0E0
            if _cs == "C":
                _mult = _scale
            for _i in range(_natoms):
                try:
                    _line = _f.readline().strip()
                    _pos.append([float(_x) for _x in _line.split()[:3]])
                except ValueError:
                    _f.close()
                    raise poscarError("Bad internal coordinates at atom line {}".format(_i+1))
            _pos = np.array(_pos, dtype=cls._dtype) * _mult
            _f.close()
            return cls(_cell, _atoms, _pos, unit="ang", coordSys=_cs, fSelectDyn=_fSelectDyn)

    @classmethod
    def create_from_lattice(cls, latt):
        '''Create POSCAR from ``lattice`` instance ``latt``.
        '''
        try:
            assert isinstance(latt, lattice)
        except AssertionError:
            raise poscarError("the input is not a lattice instance")
        __kw = latt.get_kwargs()
        return cls(*latt.get_latt(), **__kw)
    