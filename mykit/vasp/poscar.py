# -*- coding: utf-8 -*-
'''Define the class for manipulating POSCAR, the VASP lattice input
'''
import logging
import os
import string

import numpy as np

from mykit.core.cell import Cell, atoms_from_sym_nat, sym_nat_from_atoms
from mykit.core.utils import trim_after


class PoscarError(Exception):
    pass


class Poscar(Cell):
    '''The class to manipulate POSCAR, the VASP lattice input file.
    '''

    _error = PoscarError

    def __init__(self, cell, atoms, pos, **kwargs):
        super(Poscar, self).__init__(cell, atoms, pos, **kwargs)

    # ? Rewrite atomType and typeIndex properties, as the element symbol can appear twice in POSCAR
    # ? But it may be good to let user deal with it, such as using "Fe1", "Fe2" to distinguish.
    # ? In this case, it should be careful to set POTCAR when recognizing atomic information in POSCAR

    def __str__(self):
        return self.__str()

    def __str(self, scale=1.0):
        '''return the POSCAR lines in one string
        '''
        _ret = []
        # convert to ang, as vasp use ang only
        _uWas = self.unit
        self.unit = "ang"
        _syms, _nats = sym_nat_from_atoms(self.atoms)
        _ret.append(self.comment)
        _ret.append("{:8.6f}".format(scale))
        for i in range(3):
            _ret.append("  %12.8f  %12.8f  %12.8f"
                        % (self.latt[i, 0], self.latt[i, 1], self.latt[i, 2]))
        if not _syms[0].startswith("Unk"):
            _ret.append(' '.join(_syms))
        _ret.append(' '.join([str(x) for x in _nats]))
        if self.useSelDyn:
            _ret.append("Selective Dynamics")
        _ret.append({"D": "Direct", "C": "Cart"}[self.coordSys])
        for i in range(self.natoms):
            if self.useSelDyn:
                _dyn = self.sdFlags(ia=i)
            else:
                _dyn = []
            if not _syms[0].startswith("Unk"):
                _ainfo = ['#{}'.format(self.atoms[i])]
            else:
                _ainfo = []
            _aflag = [{True: "T", False: "F"}[_d] for _d in _dyn] + _ainfo
            _ret.append("%15.9f %15.9f %15.9f " % (
                self.pos[i, 0], self.pos[i, 1], self.pos[i, 2]) + ' '.join(_aflag))
        # convert back to the original length unit
        self.unit = _uWas
        return '\n'.join(_ret)

    def write(self, pathPoscar='POSCAR', scale=1.0, backup=False, suffix="_bak"):
        '''Write POSCAR to path 
        '''
        # TODO test it
        _name = pathPoscar
        try:
            assert not os.path.isdir(_name)
        except AssertionError:
            raise self._error("The path to write POSCAR is a directory.")
        if os.path.isfile(_name) and backup:
            _bakname = _name + suffix.strip()
            os.rename(_name, _bakname)
        with open(_name, 'w') as f:
            print(self.__str(scale=scale), file=f)

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
        _fAtomsAtHead = False
        try:
            _f = open(pathPoscar, 'r')
        except FileNotFoundError as _err:
            raise cls._error("Fail to open file: {}".format(pathPoscar))
        else:
            _symTypes = None
            # line 1: comment on system
            _comment = _f.readline().strip()
            # line 2: scale
            _scale = float(_f.readline().strip())
            # line 3-5: lattice vector
            _latt = [_f.readline().split() for i in range(3)]
            try:
                _latt = np.array(_latt, dtype=cls._dtype) * _scale
            except ValueError:
                _f.close()
                raise cls._error("Bad lattice vector: {}".format(pathPoscar))
            # Next 2 or 1 line(s), depend on whether element symbols are typed or not
            _line = _f.readline().strip()
            if _line[0] in string.ascii_letters:
                _fAtomsAtHead = True
                _symTypes = _line.split()
                _line = _f.readline().strip()
            if _line[0] in string.digits[1:]:
                _natomsType = [int(_x) for _x in _line.split()]
                if _symTypes is None:
                    # cls.print_cm_warn("No atom information in POSCAR: {}".format(pathPoscar))
                    _symTypes = ["Unk{}".format(i)
                                 for i, _x in enumerate(_natomsType)]
            else:
                _f.close()
                raise cls._error(
                    "Bad POSCAR atomic format: {}".format(pathPoscar))
            try:
                assert len(_symTypes) == len(_natomsType)
            except AssertionError:
                _f.close()
                raise cls._error(
                    "Inconsistent input of symbol and numbers of atom: {}".format(pathPoscar))
            _atoms = atoms_from_sym_nat(_symTypes, _natomsType)
            _natoms = sum(_natomsType)
            # Next 2 or 1 line(s), depend on whether 'selective dynamics line' is typed
            _line = _f.readline().strip()
            if _line[0].upper() == "S":
                # __flagSelDyn = True
                _line = _f.readline().strip()
            if _line[0].upper() in ["C", "K", "D"]:
                _cs = {"C": "C", "K": "C", "D": "D"}[_line[0].upper()]
            else:
                raise cls._error(
                    "Bad coordinate system: {}".format(pathPoscar))
            # Next _natoms lines: read atomic position and selective dynamics flag
            _pos = []
            _mult = 1.0E0
            _atomsFromPosLine = []
            if _cs == "C":
                _mult = _scale
            for _i in range(_natoms):
                # read the positions
                try:
                    _line = _f.readline().strip()
                    _words = trim_after(_line, r'[\#\!]').split()
                    _pos.append([float(_x) for _x in _words[:3]])
                except (ValueError, IndexError):
                    _f.close()
                    raise cls._error(
                        "Bad internal coordinates at atom line {}: {}".format(_i+1, pathPoscar))
                # read possible selective dynamic flags, and atom type
                if len(_words) == 3:
                    pass
                # add possible atomic info for ATAT like POSCAR
                elif len(_words) == 4:
                    try:
                        assert _words[-1] not in __fixDict
                    except AssertionError:
                        _f.close()
                        raise cls._error(
                            "Bad selective dynamics flag at atom line {}: {}".format(_i+1, pathPoscar))
                    _atomsFromPosLine.append(_words[-1])
                elif len(_words) in [6, 7]:
                    __fixFlag = [__fixDict.get(_words[i]) for i in range(3, 6)]
                    if None in __fixFlag:
                        # raise IndexError
                        _f.close()
                        raise cls._error(
                            "Bad selective dynamics flag at atom line {}: {}".format(_i+1, pathPoscar))
                    elif __fixFlag == [True, True, True]:
                        pass
                    else:
                        __fix.update({_i: __fixFlag})
                    if len(_words) == 7:
                        _atomsFromPosLine.append(_words[-1])
                else:
                    _f.close()
                    raise cls._error(
                        "Bad column numbers at atom line {}: {}".format(_i+1, pathPoscar))
            _pos = np.array(_pos, dtype=cls._dtype) * _mult
            _f.close()
            # print(_fAtomsAtHead, _atomsFromPosLine)
            if not _fAtomsAtHead and _atomsFromPosLine != []:
                _atoms = _atomsFromPosLine
            return cls(_latt, _atoms, _pos, unit="ang", coordSys=_cs,
                       allRelax=True, selectDyn=__fix, comment=_comment)

    @classmethod
    def create_from_cell(cls, cell):
        '''Create POSCAR from ``Cell`` instance ``cell``.
        '''
        try:
            assert isinstance(cell, Cell)
        except AssertionError:
            raise cls._error("the input is not a lattice instance")
        __kw = cell.get_kwargs()
        return cls(*cell.get_cell(), **__kw)


# ===============================

# def read_atom_info(poscar='POSCAR'):
#     atom_type_list = []
#     natoms_list = []
#     with open(poscar,'r') as h_poscar:
#         for i in range(8):
#             line = h_poscar.readline()
#             if i == 5:
#                 atom_type_list = line.split()
#             if i == 6:
#                 natoms_list = [int(x) for x in line.split()]

#     return atom_type_list, natoms_list
