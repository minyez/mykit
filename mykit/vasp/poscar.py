# -*- coding: utf-8 -*-
'''Define the class for manipulating POSCAR, the VASP lattice input
'''
import numpy as np
import logging
import string
from mykit.core.lattice import lattice, atoms_from_sym_nat, sym_nat_from_atoms

class poscarError(Exception):
    pass

class poscar(lattice):
    '''The class to manipulate POSCAR, the VASP lattice input file.
    '''

    def __init__(self, cell, atoms, pos, **kwargs):
        super(poscar, self).__init__(cell, atoms, pos, **kwargs)

    # ? Rewrite atomType and typeIndex properties, as the element symbol can appear twice in POSCAR
    # ? Maybe good to let user deal with it, such as using "Fe1", "Fe2" to distinguish
    # ? In this case, should be careful to set POTCAR when recognizing atomic information in POSCAR

    def write(self, toPOSCAR='POSCAR'):
        '''Write POSCAR to ``toPOSCAR``
        '''
        # TODO test it
        with open(toPOSCAR, 'w') as f:
            _cell, _atoms, _pos = self.get_latt()
            _syms, _nats = sym_nat_from_atoms(_atoms)
            print(self.comment, file=f)
            print(1.00000, file=f)
            for i in range(3):
                print(_cell[i,:], file=f)
            print(_syms, file=f)
            print(_nats, file=f)
            if self.useSelDyn:
                print("Selective Dynamics", file=f)
            print({"D": "Direct", "C": "Cart"}[self.coordSys], file=f)
            for i in range(self.natoms):
                if self.useSelDyn:
                    _dyn = self.sdFlags(ia=i)
                else:
                    _dyn = []
                _dynFlag = [{True:"T", False:"F"}[_d] for _d in _dyn] + ['#{}'.format(_atoms[i])]
                print(self[i], _dynFlag, file=f)
    
    @classmethod
    def read_from_file(cls, poscarPath):
        '''Read from an existing POSCAR file ``poscarPath``

        Args :
            poscarPath (str) : path to the file to read as POSCAR
        '''
        # TODO may be decomposed to support reading structures stored in XDATCAR and OUTCAR?
        __flagSelDyn = False
        __fix = {}
        __fixDict = {'T': True, 'F': False}
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
                raise poscarError("Bad lattice vector: {}".format(poscarPath))
            # Next 2 or 1 line(s), depend on whether element symbols are typed or not
            _line = _f.readline().strip()
            if _line[0] in string.ascii_letters:
                _symTypes = _line.split()
                _line = _f.readline()
            if _line[0] in string.digits[1:]:
                _natomsType = [int(_x) for _x in _line.split()]
                if _symTypes is None:
                    cls.print_cm_warn("No atom information in POSCAR: {}".format(poscarPath))
                    _symTypes = [string.ascii_lowercase[_i] for _i,_x in enumerate(_natomsType)]
            else:
                _f.close()
                raise poscarError("Bad POSCAR format: {}".format(poscarPath))
            try:
                assert len(_symTypes) == len(_natomsType)
            except AssertionError:
                _f.close()
                raise poscarError("Inconsistent input of symbol and numbers of atom: {}".format(poscarPath))
            _atoms = atoms_from_sym_nat(_symTypes, _natomsType)
            _natoms = sum(_natomsType)
            # Next 2 or 1 line(s), depend on whether 'selective dynamics line' is typed
            _line = _f.readline().strip()
            if _line[0].upper() not in ["C", "D", "K"]:
                __flagSelDyn = True
                _line = _f.readline().strip()
            if _line[0].upper() in ["C", "K"]:
                _cs = "C"
            else:
                _cs = "D"
            # Next _natoms lines: read atomic position and selective dynamics flag
            _pos = []
            _mult = 1.0E0
            if _cs == "C":
                _mult = _scale
            for _i in range(_natoms):
                try:
                    _line = _f.readline().strip()
                    _words = _line.split()
                    # remove comment part
                    for _j, _s in enumerate(_words):
                        if _s.startswith("#"):
                            _words = _words[0:_j+1]
                            break
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
                    raise poscarError("Bad internal coordinates at atom line {}: {}".format(_i+1, poscarPath))
                except IndexError:
                    _f.close()
                    raise poscarError("Bad selective dynamics flag at atom line {}: {}".format(_i+1, poscarPath))
            _pos = np.array(_pos, dtype=cls._dtype) * _mult
            _f.close()
            return cls(_cell, _atoms, _pos, unit="ang", coordSys=_cs, allRelax=True, selectDyn=__fix)

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

    # TODO inherit factory methods from ``lattice`` for space group and Bravis cell creation
    