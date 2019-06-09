# -*- coding: utf-8 -*-
# pylint: disable=bad-whitespace
'''Module that defines classes for crystal cell manipulation and symmtetry operation

The ``cell`` class and its subclasses accept the following kwargs when being instantialized:

    - coordSys (str): Coordinate system for the internal positions,
      either "D" (Direct, default) or "C" (Cartesian)
    - allRelax (bool) : default selective dynamics option for atoms.
      Set True (default) to allow all DOFs to relax
    - selectDyn (dict) : a dictionary with key-value pair as ``int: [bool, bool, bool]``, 
      which controls the selective dynamic option for atom with the particular index 
      (starting from 0). Default is an empty ``dict``
    - comment (str): message about the cell
    - reference (str): the reference where the lattice structure is derived.

When other keyword are parsed, they will be filtered out and no exception will be raised
'''
from collections import OrderedDict
from numbers import Real

import numpy as np

from mykit.core.constants import PI
from mykit.core.log import Verbose
from mykit.core.numeric import Prec
from mykit.core.unit import LengthUnit
from mykit.core.utils import (Cif, get_latt_consts_from_latt_vecs,
                              get_str_indices)


# ==================== classes ====================
class CellError(Exception):
    '''Exception in cell module
    '''
    pass


class Cell(Prec, Verbose, LengthUnit):
    '''Cell structure class

    Args:
        latt (array-like) : The lattice vectors
        atoms (list of str) : The list of strings of type for each atom 
        corresponding to the member in pos
        pos (array-like) : The internal coordinates of atoms
        unit (str): the unit, in lower case, either "ang" (default) or "au".

    Note:
        see ``cell`` module docstring for acceptable kwargs for ``Cell`` and its subclasses

    Examples:
    >>> Cell([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], ["C"], [[0.0, 0.0, 0.0]])
    <mykit.core.cell.Cell>
    '''

    _error = CellError

    def __init__(self, latt, atoms, pos, unit='ang', **kwargs):

        self.comment = 'Default Cell class'
        self.__reference = ''
        self.__allRelax = True
        self.__selectDyn = {}
        self.__coordSys = 'D'

        try:
            self.__latt = np.array(latt, dtype=self._dtype)
            self.__pos = np.array(pos, dtype=self._dtype)
        except ValueError as _err:
            raise self._error(
                "Fail to create cell and pos array. Please check.")
        LengthUnit.__init__(self, lunit=unit)
        self.__atoms = [_a.capitalize() for _a in atoms]
        self.__parse_cellkw(**kwargs)
        # check input consistency
        self.__check_consistency()
        # move all atoms into the lattice (0,0,0)
        # self.move_atoms_to_first_lattice()
        # sanitize atoms arrangement
        self.__sanitize_atoms()

    def __len__(self):
        return len(self.__atoms)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.__pos[index, :]
        if isinstance(index, str):
            return self.get_sym_index(index)
        raise self._error("atom index or symbol {} not found".format(index))

    def __str__(self):
        return "{}\nLattice:\n{}\nAtoms: {}\nPositions:\n{}\nUnit: {}\nCoordinates in {}".format(
            self.comment, self.latt, self.atoms, self.pos, self.unit, self.coordSys)

    def __repr__(self):
        return self.__str__()

    def __parse_cellkw(self, **kwargs):
        if 'coordSys' in kwargs:
            self.__coordSys = kwargs['coordSys'].upper()
        if "allRelax" in kwargs:
            self.__allRelax = kwargs["allRelax"]
        if "selectDyn" in kwargs:
            self.__selectDyn = kwargs["selectDyn"]
        if "reference" in kwargs:
            self.__reference = "{}".format(kwargs["reference"])
        if "comment" in kwargs:
            self.comment = "{}".format(kwargs["comment"])

    def get_cell(self):
        '''Purge out the cell, atoms and pos.

        ``cell``, ``atoms`` and ``pos`` are the minimal for constructing ``Cell``
        and its subclasses.
        They can also be used to build sysmetry operations with spglib utilities.
        '''
        return self.__latt, self.__atoms, self.__pos

    def get_kwargs(self):
        '''return all kwargs useful to constract program-dependent cell input from ``Cell`` instance

        Returns :
            dictionary that can be parsed to ``create_from_cell`` class method.
        '''
        _d = {
            "unit": self._lunit,
            "coordSys": self.__coordSys,
            "comment": self.comment,
            "reference": self.__reference,
            "allRelax": self.__allRelax,
            "selectDyn": self.__selectDyn
        }
        return _d

    def get_reference(self):
        '''Return the reference of the structure
        '''
        return self.__reference

    def __check_consistency(self):
        try:
            assert self.__coordSys in ["C", "D"]
            assert np.shape(self.__latt) == (3, 3)
            assert self.natoms > 0
            assert np.shape(self.__pos) == (self.natoms, 3)
        except AssertionError:
            raise self._error("Invalid cell setup")
        # ? switch automatically, or let the user deal with it
        try:
            assert self.vol > 0
        except AssertionError:
            raise self._error(
                "Left-handed system found (vol<0). Switch two vector.")

    def _switch_two_atom_index(self, iat1, iat2):
        '''switch the index of atoms with index iat1 and iat2

        Except ``__pos`` and ``__atoms``,
        this method should also deals possible switch in other positional
        attributes, e.g.

        - ``selectDyn`` (DONE)

        Note that this method is mainly for sorting use, and does NOT change
        the geometry of the cell at all.
        '''
        try:
            assert iat1 in range(self.natoms)
            assert iat2 in range(self.natoms)
            assert iat1 != iat2
        except AssertionError:
            raise self._error(
                "Fail to switch two atoms with indices {} and {}".format(iat1, iat2))

        self.__pos[[iat1, iat2]] = self.__pos[[iat2, iat1]]
        self.__atoms[iat1], self.__atoms[iat2] = self.__atoms[iat2], self.__atoms[iat1]

        _sfd1 = self.__selectDyn.pop(iat1, [])
        _sfd2 = self.__selectDyn.pop(iat2, [])
        if _sfd1 != []:
            self.__selectDyn.update({iat2: _sfd1})
        if _sfd2 != []:
            self.__selectDyn.update({iat1: _sfd2})

    def move_atoms_to_first_lattice(self):
        '''Move all atoms into the lattice (0,0,0).

        For Cartisian system, the move is achieved by first
        converting to and then back from direct system.
        By this process, each component of the coordinate (in direct system)
        belongs to [0,1)
        '''
        if self.coordSys == "D":
            self.__pos = self.__pos - np.floor(self.__pos)
        elif self.coordSys == "C":
            self.coordSys = "D"
            self.__pos = self.__pos - np.floor(self.__pos)
            self.coordSys = "C"

    def get_sym_index(self, csymbol):
        '''Get the indices of atoms with element symbol ``csymbol``

        Note that this is equivalent to ``cell[csymbol]``, given cell an instance of ``Cell``.

        Args:
            csymbol (str) : chemical-symbol-like identifier
        '''
        assert isinstance(csymbol, str)
        return get_str_indices(self.atoms, csymbol.capitalize())

    # * Sorting method
    def __bubble_sort_atoms(self, key, indices, reverse=False):
        '''sort atoms with bubble sort under various scenarios

        The smaller value will appear earlier, if ``reverse`` is left
        as False.
        In both cases, when two same values are compared,
        current bubble will just break.

        Args:
            key (natom-member list): the key value to be sorted
            indices (iterable): the indices of the atoms to be sorted
            reverse (bool): if set True, larger value appears earlier
        '''
        _depth = 1
        # self.print_log("Bubble sort with key:", key, ", indices:", indices, level=3, depth=_depth)
        __ind = list(indices)
        __key = [key[_i] for _i in __ind]
        _n = len(__ind)
        __sorted = True
        for _i in range(_n - 1):
            _li = _i
            _ri = _i+1
            # self.print_log("Check index: ", _i, level=3, depth=_depth+2)
            __dict = {True: __key[_li] > __key[_ri],
                      False: __key[_li] < __key[_ri]}
            if not __dict[reverse]:
                __sorted = False
                break
        if not __sorted:
            for _i in range(1, _n):
                _j = _i
                # self.print_log("Sorting size {}".format(_i+1), level=3, depth=_depth+1)
                while _j > 0:
                    _li = _j-1
                    _ri = _j
                    __dict = {True: __key[_ri] > __key[_li],
                              False: __key[_ri] < __key[_li]}
                    if __dict[reverse]:
                        self._switch_two_atom_index(__ind[_li], __ind[_ri])
                        __key[_li], __key[_ri] = __key[_ri], __key[_li]
                        _j -= 1
                    else:
                        break
                # self.print_log("Sorting size {} done".format(_i+1), level=3, depth=_depth+1)
        # self.print_log("Bubble sort done", level=3, depth=_depth)

    def __sanitize_atoms(self):
        '''Sanitize the atoms arrangement after initialization.

        It mainly deals with arbitrary input of ``atoms`` when initialized.
        '''
        self.__bubble_sort_atoms(self.typeIndex, range(self.natoms))

    def sort_pos(self, axis=3, reverse=False):
        '''Sort the atoms by its coordinate along axis.

        The ``atoms`` list will not change by sorting, i.e. the sorting is performed
        within each atomic group.
        If ``reverse`` is False, atom with higher coordinate in lattice (0,0,0)
        will appear earlier, otherwise later.

        This behavior is opposite to ``sort`` functions, in spirit of that surfaces are often
        placed at high along one axis, and sorting in descending order makes the surface
        appear first and easy to modify.

        Args :
            axis (1,2,3)
            reverse (bool)
        '''
        try:
            assert axis in range(1, 4)
            assert isinstance(reverse, bool)
        except AssertionError:
            raise self._error()
        __sortKeys = self.pos[:, axis-1]
        # self.print_log("__sortKeys in sort_pos", __sortKeys, level=3)
        for _at in self.atomTypes:
            __ind = self.get_sym_index(_at)
            self.__bubble_sort_atoms(__sortKeys, __ind, reverse=not reverse)

    # * Cell manipulation
    def scale(self, scale):
        '''Scale the lattice, i.e. increase the lattice by ``scale`` time
        '''
        try:
            assert isinstance(scale, Real)
            assert scale > 0.0
        except AssertionError:
            raise self._error("scale must be positive real")
        self.__latt = self.__latt * scale
        if self.__coordSys == "C":
            self.__pos = self.__pos * scale

    def add_atom(self, atom, coord, sdFlag=None):
        '''Add an atom with coordinate and selective dynamic flags

        Args:
            atom (str): the chemical symbol of the atom to add
            coord (array-like): the coordinate of atom in ``Cell`` coordinate system
            sdFlag (list of 3 bools): 
        '''
        try:
            assert isinstance(atom, str)
        except:
            raise self._error(
                "atom should be string, received {}".format(type(atom)))
        try:
            newPos = np.vstack([self.__pos, coord])
        except ValueError:
            raise self._error("Invalid coordinate: {}".format(coord))
        if sdFlag is not None:
            self.__set_sdFlags({self.natoms: sdFlag})
        self.__pos = newPos
        self.__atoms.append(atom)
        self.move_atoms_to_first_lattice()
        self.__sanitize_atoms()

    # TODO move atom
    def __move(self, ia):
        raise NotImplementedError

    def __move_all(self, shift):
        '''Move all atoms by a shift
        '''
        assert np.shape(shift) == (3,)
        np.add(self.__pos, shift, out=self.__pos)

    @property
    def center(self):
        '''Calculate the center of all atoms in the cell
        '''
        assert self.coordSys == "D"
        _posSum = np.zeros(3, dtype=self._dtype)
        _n = 0
        for i in range(self.natoms):
            _dupcs, _dn = periodic_duplicates_in_cell(self.__pos[i, :])
            for _dupc in _dupcs:
                np.add(_posSum, _dupc/float(_dn), _posSum)
        # _posSum = np.sum(self.__pos, axis=0)
        # check periodic duplicate, by recognizing number of zeros in pos.
        # _dup = 3 - np.count_nonzero(self.__pos, axis=1)
        # _dup = np.power(2, _dup)
        # _n = np.sum(_dup)
        # _posSum = np.sum(self.__pos * _dup[:, None], axis=0)
        return _posSum / self.natoms

    def centering(self, axis=0):
        '''Centering the atoms along axes. Mainly use for slab model.

        TODO:
            For now not work when there is atom at origin along the axis

        Args:
            axis (int or iterable of int) : the axes along which the atoms will be centered.
        '''
        _aList = axis_list(axis)
        _wasCart = self.coordSys == "C"
        if _wasCart:
            self.coordSys = "D"
        # get the geometric center of all atoms
        _center = self.center
        _shift = np.array([0.5, 0.5, 0.5], dtype=self._dtype) - _center

        for i in range(3):
            ia = i + 1
            if ia not in _aList:
                _shift[i] = 0.0
        self.__move_all(_shift)
        #     if self.check_vacuum_pos(zdirt):
        #         self.__print("  - Vacuum in the middle detected. Not supported currently. Pass.")
        #         continue
        #     else:
        #         surf_atom = [self.check_extreme_atom(0.0,zdirt,False,1.0),self.check_extreme_atom(0.0,zdirt,True,1.0)]
        #         # debug
        #         # print surf_atom
        #         shift = 0.5 - sum([self.innerpos[i-1][iz] for i in surf_atom])/2.0
        #         self.action_shift(shift,zdirt)
        # self.__print(" Complete centering.")
        if _wasCart:
            self.coordSys = "C"

    @property
    def a(self):
        '''Lattice vectors
        '''
        return self.__latt

    @property
    def alen(self):
        '''Length of lattice vectors
        '''
        return np.array([np.linalg.norm(x) for x in self.__latt], dtype=self._dtype)

    @property
    def lattConsts(self):
        '''Lattice constant of the cell, i.e., a, b, c, alpha, beta, gamma (in degree)
        '''
        return get_latt_consts_from_latt_vecs(self.__latt)

    @property
    def latt(self):
        '''Lattice vectors
        '''
        return self.__latt

    @property
    def atoms(self):
        '''list.'''
        return self.__atoms

    @property
    def pos(self):
        '''array.'''
        return self.__pos

    @property
    def unit(self):
        '''str.'''
        return self._lunit

    @unit.setter
    def unit(self, u):
        coef = self._get_lunit_conversion(u)
        if coef != 1:
            if self.__coordSys == "C":
                self.__pos = self.__pos * coef
            self.__latt = self.__latt * coef
            self._lunit = u.lower()

    @property
    def coordSys(self):
        '''coordinate system
        '''
        return self.__coordSys

    @coordSys.setter
    def coordSys(self, s):
        _s = s.upper()
        if _s != self.__coordSys:
            _convDict = {"C": self.__latt, "D": np.linalg.inv(self.__latt)}
            _conv = _convDict.get(_s)
            if _conv is not None:
                self.__pos = np.matmul(self.__pos, _conv)
                self.__coordSys = _s
            else:
                info = "Only support \"D\" direct or fractional and \"C\" Cartisian coordinate."
                raise self._error(info)

    @property
    def atomTypes(self):
        '''All atom types in the cell
        '''
        _d = OrderedDict.fromkeys(self.__atoms)
        return list(_d.keys())

    @property
    def typeMapping(self):
        '''Map index (int) to atom type (str)
        '''
        _ats = self.atomTypes
        _dict = {}
        for i, _at in enumerate(_ats):
            _dict.update({i: _at})
        return _dict

    @property
    def typeIndex(self):
        '''Indices of atom type of all atoms
        '''
        _ats = self.atomTypes
        _dict = {}
        for i, _at in enumerate(_ats):
            _dict.update({_at: i})
        return [_dict[_a] for _a in self.__atoms]

    @property
    def vol(self):
        '''Volume of the cell
        '''
        return np.linalg.det(self.__latt)

    @property
    def natoms(self):
        '''Int. Total number of atoms
        '''
        return len(self.__atoms)

    @property
    def useSelDyn(self):
        if self.__allRelax and not bool(self.__selectDyn):
            return False
        return True

    @property
    def recpLattIn2Pi(self):
        '''Reciprocal lattice vectors in 2Pi unit^-1
        '''
        b = []
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            b.append(np.cross(self.latt[j, :], self.latt[k, :]))
        return np.array(b, dtype=self._dtype) / self.vol

    @property
    def bIn2Pi(self):
        '''Alias to ``recpLattIn2Pi``
        '''
        return self.recpLattIn2Pi

    @property
    def recpLatt(self):
        '''Reciprocal lattice vectors in unit^-1
        '''
        return self.bIn2Pi * 2.0E0 * PI

    @property
    def b(self):
        '''Alias of ``recpLatt``
        '''
        return self.recpLatt

    @property
    def blen(self):
        '''Length of reciprocal lattice vector in unit^-1
        '''
        return np.array([np.linalg.norm(x) for x in self.b], dtype=self._dtype)

    # * selective dynamics related
    def fix_all(self):
        '''Fix all atoms.
        '''
        self.__selectDyn = {}
        self.__allRelax = False

    def relax_all(self):
        '''Relax all atoms.
        '''
        self.__selectDyn = {}
        self.__allRelax = True

    def set_fix(self, *iats, axis=0):
        '''Fix the atoms with index in iats

        Args:
            iats (list of int): the indices of atoms to fix
            axis (int or list): the axes along which the position of atom is fixed
                It can be 0|1|2|3, or a list with all its members 1|2|3
        '''
        if len(iats) != 0:
            _new = {}
            for _ia in iats:
                if _ia in range(self.natoms):
                    _new.update(
                        {_ia: select_dyn_flag_from_axis(axis, relax=False)})
            self.__set_sdFlags(_new)

    def set_relax(self, *iats, axis=0):
        '''Relax the atoms with index in iats

        Args:
            iats (list of int): the indices of atoms to relax
            axis (int or list): the axes along which the position of atom is relaxed
                It can be 0|1|2|3, or a list with all its members 1|2|3
        '''
        if len(iats) != 0:
            _new = {}
            for _ia in iats:
                if _ia in range(self.natoms):
                    _new.update(
                        {_ia: select_dyn_flag_from_axis(axis, relax=True)})
            self.__set_sdFlags(_new)

    def relax_from_top(self, n, axis=3):
        '''Set all atoms fixed, and relax the n atoms from top along axis
        '''
        pass

    def fix_from_center(self, n, axis=3):
        '''Set all atoms relaxed, and fix the n atoms from the middle along axis
        '''
        pass

    def __set_sdFlags(self, selectDyn):
        try:
            assert isinstance(selectDyn, dict)
        except AssertionError:
            raise self._error("need dictionary to set selective dynamics")
        for _k, flag in selectDyn.items():
            try:
                assert isinstance(flag, list)
                assert len(flag) == 3
                assert all([isinstance(_x, bool) for _x in flag])
            except AssertionError:
                raise self._error("Bad flag for selective dynamics")
        self.__selectDyn.update(selectDyn)

    def sdFlags(self, ia=-1):
        '''Return the selective dynamic flag (bool) of atom

        Args:
            ia (int) : index of atom

        Returns:
            3-member list, if ia is in range(natoms),
            otherwise natoms-member list, each member a 3-member list
            as the flag for that atom
        '''
        # self.print_log("Global selective dynamics flag: {}".format(self.__allRelax), level=3)
        if ia in self.__selectDyn:
            # self.print_log("Found custom flag in __selectDyn for atom {}".format(ia), depth=1, level=3)
            _flag = self.__selectDyn[ia]
        elif ia in range(self.natoms):
            # self.print_log("Use global flag for atom {}".format(ia), level=3, depth=1)
            _flag = [self.__allRelax, ] * 3
        else:
            _flag = [[self.__allRelax, ]*3 for _i in range(self.natoms)]
            for i in self.__selectDyn:
                _flag[i] = self.__selectDyn[i]
        return _flag

    def get_spglib_input(self):
        '''Return the input necessary for spglib to get symmetry

        Returns:
            cell (3,3), pos (n,3), index of atom type (n), with n = self.natoms
        '''
        return self.latt, self.pos, self.typeIndex

    # * Factory methods
    @classmethod
    def read_from_json(cls, pathJson):
        '''Initialize a ``Cell`` instance from a JSON file

        If "factory" key does not exist, it will search for the postional arguments,
        i.e. "latt", "atoms" and "pos" keys. Raise when any of them does not exist.

        Args:
            pathJson (str): the path of JSON file
        '''
        import json
        import os
        if pathJson is None or not os.path.isfile(pathJson):
            raise cls._error("JSON file not found: {}".format(pathJson))
        with open(pathJson, 'r') as h:
            try:
                js = json.load(h)
            except json.JSONDecodeError:
                raise cls._error(
                    "invalid JSON file for cell: {}".format(pathJson))
        pargs = []
        factoryDict = {
            "bravais_oP": (cls.bravais_oP, ("atom", "a", "b", "c")),
            "bravais_oI": (cls.bravais_oI, ("atom", "a", "b", "c")),
            "bravais_oF": (cls.bravais_oF, ("atom", "a", "b", "c")),
            "bravais_cP": (cls.bravais_cP, ("atom", "a")),
            "bravais_cI": (cls.bravais_cI, ("atom", "a")),
            "bravais_cF": (cls.bravais_cF, ("atom", "a")),
            "perovskite": (cls.perovskite, ("atom1", "atom2", "atom3", "a")),
            "zincblende": (cls.zincblende, ("atom1", "atom2", "a")),
            "diamond": (cls.diamond, ("atom", "a")),
            "wurtzite": (cls.wurtzite, ("atom1", "atom2", "a")),
            "rutile": (cls.rutile, ("atom1", "atom2", "a", "c", "u")),
            "anatase": (cls.anatase, ("atom1", "atom2", "a", "c", "u")),
            "pyrite": (cls.pyrite, ("atom1", "atom2", "a", "u")),
            "marcasite": (cls.marcasite, ("atom1", "atom2", "a", "b", "c", "v", "w")),
        }
        # found factory key
        if "factory" in js:
            fac = js["factory"]
            # pop out latt, atoms and pos for safety
            for arg in ["latt", "atoms", "pos"]:
                js.pop(arg, None)
            if fac in factoryDict:
                # get positional argument
                # print(fac)
                try:
                    m, reqPa = factoryDict[fac]
                    for x in reqPa:
                        pargs.append(js.pop(x))
                    return m(*pargs, **js)
                except KeyError:
                    raise cls._error(
                        "Required key not found in JSON: {}".format(x))
            else:
                raise cls._error("Factory method unavailable: {}".format(fac))

        for _i, arg in enumerate(["latt", "atoms", "pos"]):
            v = js.pop(arg, None)
            if v is None:
                raise cls._error(
                    "invalid JSON file for cell: {}. No {}".format(pathJson, arg))
            pargs.append(v)
        return cls(*pargs, **js)

    @classmethod
    def read_from_cif(cls, pathCif):
        '''Read from Cif file and return a instance by use of PyCIFRW
        '''
        cif = Cif(pathCif)
        kw = {"coordSys": "D", "reference": cif.get_reference_str(), }
        # use chemical name as comment
        kw['comment'] = ', '.join(cif.get_chemical_name()) + ' type'
        latt = cif.get_lattice_vectors()
        atoms, pos = cif.get_all_atoms()
        return cls(latt, atoms, pos, **kw)

    @classmethod
    def create_from_cell(cls, cell):
        '''Create an ``Cell`` object from another ``Cell`` instance ``cell``.

        This is for use of transformation between cells described 
        by different file formats.

        Args:
            cell : object of ``Cell`` or its subclasses

        Returns:
            a ``Cell`` object
        '''
        try:
            assert isinstance(cell, Cell)
        except AssertionError:
            raise cls._error(
                "the input is not an object of Cell or its subclasses")
        kw = cell.get_kwargs()
        return cls(*cell.get_cell(), **kw)

    @classmethod
    def _bravais_o(cls, kind, atom, a, b, c, **kwargs):
        assert kind in ["P", "I", "F"]
        _latt = [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
        if kind == "P":
            _atoms = [atom, ]
            _pos = [[0.0, 0.0, 0.0]]
        if kind == "I":
            _atoms = [atom, ]*2
            _pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        if kind == "F":
            _atoms = [atom, ]*4
            _pos = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        kwargs.pop("coordSys", None)
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def bravais_oP(cls, atom, a=1.0, b=2.0, c=3.0, **kwargs):
        '''Generate a simple orthorhombic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            a,b,c (float) : the lattice constants (a,b,c)
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        if "comment" not in kwargs:
            kwargs.update(
                {"comment": "Simple orthorhombic lattice {}".format(atom)})
        return cls._bravais_o("P", atom, a, b, c, **kwargs)

    @classmethod
    def bravais_oI(cls, atom, a=1.0, b=2.0, c=3.0, **kwargs):
        '''Generate a body-centered orthorhombic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            a,b,c (float) : the lattice constants (a,b,c)
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        if "comment" not in kwargs:
            kwargs.update(
                {"comment": "Body-centered orthorhombic lattice {}".format(atom)})
        return cls._bravais_o("I", atom, a, b, c, **kwargs)

    @classmethod
    def bravais_oF(cls, atom, a=1.0, b=2.0, c=3.0, **kwargs):
        '''Generate a face-centered orthorhombic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            a,b,c (float) : the lattice constants (a,b,c)
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        if "comment" not in kwargs:
            kwargs.update(
                {"comment": "Face-centered orthorhombic lattice {}".format(atom)})
        return cls._bravais_o("F", atom, a, b, c, **kwargs)

    @classmethod
    def bravais_cP(cls, atom, a=1.0, **kwargs):
        '''Generate a simple cubic Bravais lattice, space group 221

        Args:
            atom (str) : the chemical symbol of atom
            a (float) : the lattice constant (a)
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _a]]
        _atoms = [atom]
        _pos = [[0.0, 0.0, 0.0]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "Simple cubic lattice {}".format(atom)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def bravais_cI(cls, atom, a=1.0, primitive=False, **kwargs):
        '''Generate a body-centered cubic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            a (float) : the lattice constant (a)
            primitive (bool) : if set True, the primitive cell will be generated
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        if primitive:
            _latt = [[-_a/2.0, _a/2.0, _a/2.0],
                     [_a/2.0, -_a/2.0, _a/2.0], [_a/2.0, _a/2.0, -_a/2.0]]
            _atoms = [atom]
            _pos = [[0.0, 0.0, 0.0]]
        else:
            _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _a]]
            _atoms = [atom, ]*2
            _pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "BCC {}".format(atom)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def bravais_cF(cls, atom, a=1.0, primitive=False, **kwargs):
        '''Generate a face-centered cubic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            a (float) : the lattice constant (a)
            primitive (bool) : if set True, the primitive cell will be generated
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        if primitive:
            _latt = [[0.0, _a/2.0, _a/2.0],
                     [_a/2.0, 0.0, _a/2.0], [_a/2.0, _a/2.0, 0.0]]
            _atoms = [atom]
            _pos = [[0.0, 0.0, 0.0]]
        else:
            _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _a]]
            _atoms = [atom, ]*4
            _pos = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "FCC {}".format(atom)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def perovskite(cls, atom1="Ca", atom2="Ti", atom3="O", a=1.0, **kwargs):
        '''Generate a perovskit lattice

        Args:
            atom1 (str) : the chemical symbol of atom at vertices of cubic cell
            atom2 (str) : the chemical symbol of atom at center of cubic cell
            atom3 (str) : the chemical symbol of atom at faces of cubic cell
            a (float) : the lattice constant (a)
            primitive (bool) : if set True, the primitive cell will be generated
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _a]]
        _atoms = [atom1, atom2, ] + [atom3, ]*3
        _pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [
            0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update(
                {"comment": "Perovskite {}{}{}3".format(atom1, atom2, atom3)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def zincblende(cls, atom1="Zn", atom2="O", a=1.0, primitive=False, **kwargs):
        '''Generate a zincblende lattice (space group 216)

        ``atom1`` are placed at vertex and ``atom2`` at tetrahedron interstitial

        Args:
            atom1 (str): symbol of atom at vertex
            atom2 (str): symbol of atom at tetrahedron interstitial
            a (float): the lattice constant of the conventional cell.
            primitive (bool): if set True, the primitive cell will be generated.
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        if primitive:
            _latt = [[0.0, _a/2.0, _a/2.0],
                     [_a/2.0, 0.0, _a/2.0], [_a/2.0, _a/2.0, 0.0]]
            _atoms = [atom1, atom2]
            _pos = [[0.0, 0.0, 0.0],
                    [0.25, 0.25, 0.25]]
        else:
            _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _a]]
            _atoms = [atom1, ]*4 + [atom2, ]*4
            _pos = [[0.0, 0.0, 0.0],
                    [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5],
                    [0.5, 0.5, 0.0],
                    [0.25, 0.25, 0.25],
                    [0.25, 0.75, 0.75],
                    [0.75, 0.25, 0.75],
                    [0.75, 0.75, 0.25]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "Zincblende {}{}".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def diamond(cls, atom="C", a=1.0, primitive=False, **kwargs):
        '''Generate a diamond lattice (space group 227)

        Args:
            atom (str): symbol of the atom
            a (float): the lattice constant of the conventional cell.
            primitive (bool): if set True, the primitive cell will be generated.
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        if "comment" not in kwargs:
            kwargs.update({"comment": "Diamond {}".format(atom)})
        return cls.zincblende(atom, atom, a=a, primitive=primitive, **kwargs)

    @classmethod
    def wurtzite(cls, atom1="Zn", atom2="O", a=1.0, **kwargs):
        '''Generate a wurtzite lattice (space group 186)

        Args:
            atom1 (str): symbol of atom at vertices of lattice
            atom2 (str): symbol of atom at edges of lattice
            a (float): the lattice constant of the cell.
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        _halfa = _a/2.0
        _c = _a * np.sqrt(8.0/3)
        _latt = [[_a, 0.0, 0.0],
                 [-_halfa, np.sqrt(3)*_halfa, 0.0], [0.0, 0.0, _c]]
        _atoms = [atom1, ]*2 + [atom2, ]*2
        _pos = [[0.0, 0.0, 0.0],
                [2.0/3, 1.0/3, 0.5],
                [0.0, 0.0, 2.0/3],
                [2.0/3, 1.0/3, 1.0/6]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "Wurtzite {}{}".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def rutile(cls, atom1="Ti", atom2="O", a=1.0, c=2.0, u=0.31, **kwargs):
        '''Generate a rutile lattice (space group 136)

        Args:
            atom1 (str): symbol of atom at vertex and center of lattice
            atom2 (str): symbol of atom at face of lattice
            a,c (float): the lattice constant of the cell.
            u (float): the internal coordinate
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        try:
            assert 0.0 < u < 1.0
        except AssertionError:
            raise cls._error(
                "Internal coordinate should be in (0,1), get {}".format(u))
        _a = abs(a)
        _c = abs(c)
        _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _c]]
        _atoms = [atom1, ]*2 + [atom2, ]*4
        _pos = [[0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [u, u, 0.0],
                [-u, -u, 0.0],
                [0.5-u, 0.5+u, 0.5],
                [0.5+u, 0.5-u, 0.5]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "Rutile {}{}2".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def anatase(cls, atom1="Ti", atom2="O", a=3.7845, c=9.5143, u=0.2199,
                primitive=False, **kwargs):
        '''Generate an anatase lattice (space group 141).

        Note:
            This cell is not standardized.

        Args:
            atom1 (str): symbol of atom at vertex
            atom2 (str): symbol of atom at face of lattice
            a,c (float): the lattice constant of the conventional cell.
            u (float): the internal coordinate, i.e. distance between two atoms in terms of c
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        try:
            assert 0.0 < u < 1.0
        except AssertionError:
            raise cls._error(
                "Internal coordinate should be in (0,1), get {}".format(u))
        _a = abs(a)
        _c = abs(c)
        if primitive:
            _latt = [[-_a/2, _a/2, _c/2],
                     [_a/2, -_a/2, _c/2], [_a/2, _a/2, -_c/2]]
            _atoms = [atom1, ]*2 + [atom2, ]*4
            _pos = [[0.0, 0.0, 0.0],
                    [0.75, 0.25, 0.5],
                    [0.25-u, 0.75-u, 0.5],
                    [0.25+u, 0.75+u, 0.5],
                    [0.5+u,  0.5+u, 0.0],
                    [0.5-u,  0.5-u, 0.0], ]
        else:
            _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _c]]
            _atoms = [atom1, ]*4 + [atom2, ]*8
            _pos = [[0.0, 0.0, 0.0],
                    [0.5, 0.0, 0.25],
                    [0.0, 0.5, 0.75],
                    [0.5, 0.5, 0.5],
                    [0.0, 0.0,   u],
                    [0.0, 0.0,  -u],
                    [0.5, 0.0, 0.25-u],
                    [0.5, 0.0, 0.25+u],
                    [0.0, 0.5, 0.75-u],
                    [0.0, 0.5, 0.75+u],
                    [0.5, 0.5, 0.5-u],
                    [0.5, 0.5, 0.5+u]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "Anatase {}{}2".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def pyrite(cls, atom1="Fe", atom2="S", a=5.4183, u=0.1174, **kwargs):
        '''Generate a standardized pyrite lattice (space group 205).

        Args:
            atom1 (str): symbol of atom at vertex and face-center
            atom2 (str): symbol of atom at edges
            a (float): the lattice constant of the conventional cell.
            u (float): the internal coordinate
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _a]]
        _atoms = [atom1, ]*4 + [atom2, ]*8
        _pos = [[0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
                [0.5-u,     u,    -u],
                [0.5+u,    -u,     u],
                [-u, 0.5-u,     u],
                [u, 0.5+u,    -u],
                [u,    -u, 0.5-u],
                [-u,     u, 0.5+u],
                [0.5+u, 0.5+u, 0.5+u],
                [0.5-u, 0.5-u, 0.5-u]]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "Pyrite {}{}2".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def marcasite(cls, atom1="Fe", atom2="S",
                  a=4.4450, b=5.4151, c=3.3922,
                  v=0.2066, w=0.3750, **kwargs):
        '''Generate a standardized marcasite lattice (space group 58).

        Args:
            atom1 (str): symbol of atom at vertex and body-center
            atom2 (str): symbol of the other atom
            a, b, c (float): the lattice constants of the cell.
            v, w(float): the internal coordinates
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        _b = abs(b)
        _c = abs(c)
        _latt = [[_a, 0.0, 0.0], [0.0, _b, 0.0], [0.0, 0.0, _c]]
        _atoms = [atom1, ]*2 + [atom2, ]*4
        _pos = [[0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [0.5+v, 0.5-w,    0.5],
                [0.5-v, 0.5+w,    0.5],
                [-v,    -w,    0.0],
                [v,     w,    0.0], ]
        kwargs.pop("coordSys", None)
        if "comment" not in kwargs:
            kwargs.update({"comment": "Marcasite {}{}2".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)


def atoms_from_sym_nat(syms, nats):
    '''Generate ``atom`` list for ``Cell`` initilization from list of atomic symbols 
    and number of atoms

    Args :
        syms (list of str) : atomic symbols
        nats (list of int) : number of atoms for each symbol

    Returns :
        a list of str, containing symbol of each atom in the cell

    Examples:
    >>> atoms_from_sym_nat(["C", "Al", "F"], [2, 3, 1])
    ["C", "C", "Al", "Al", "Al", "F"]
    '''
    assert len(syms) == len(nats)
    _list = []
    for _s, _n in zip(syms, nats):
        _list.extend([_s, ] * _n)
    return _list


def sym_nat_from_atoms(atoms):
    '''Generate lists of atomic symbols and number of atoms from whole atoms list

    The order of appearence of the element is conserved in the output.

    Args :
        atoms (list of str) : symbols of each atom in the cell

    Returns :
        list of str : atomic symbols
        list of int : number of atoms for each symbol

    Examples:
    >>> sym_nat_from_atoms(["C", "Al", "Al", "C", "Al", "F"])
    ["C", "Al", "F"], [2, 3, 1]
    '''
    _syms = []
    _natsDict = {}
    for _at in atoms:
        if _at in _syms:
            _natsDict[_at] += 1
        else:
            _syms.append(_at)
            _natsDict.update({_at: 1})
    return _syms, [_natsDict[_at] for _at in _syms]


def select_dyn_flag_from_axis(axis, relax=False):
    '''Generate selective dynamic flags, i.e. [bool, bool, bool]

    Args:
        relax (bool): if True, the flag for axis will be set as True.
            Otherwise False
    '''
    assert isinstance(relax, bool)
    _flag = [not relax, not relax, not relax]
    _aList = axis_list(axis)
    for _a in _aList:
        _flag[_a-1] = not _flag[_a-1]
    return _flag


def axis_list(axis):
    '''Generate axis indices from ``axis``

    Args:
        axis (int or list of int)

    Returns:
        tuple
    '''
    _aList = []
    if isinstance(axis, int):
        if axis == 0:
            _aList = [1, 2, 3]
        if axis in range(1, 4):
            _aList = [axis]
    elif isinstance(axis, (list, tuple)):
        _aSet = list(set(axis))
        for _a in _aSet:
            try:
                assert isinstance(_a, int)
            except AssertionError:
                pass
            else:
                if _a == 0:
                    _aList = [1, 2, 3]
                    break
                if _a in range(1, 4):
                    _aList.append(_a)
    return tuple(_aList)


def periodic_duplicates_in_cell(directCoord):
    '''Return the coordinates and numbers of the duplicates of an atom
    in a cell due to lattice translation symmetry

    Args:
        directCoord (array): the direct coordinate of an atom in the cell

    Note:
        The function works only when each component belongs to [0,1)

    TODO:
        Generalize this function to mirrors in n-th lattice shell

    Returns:
        tuple : the coordinates of all the atom duplicates due to transilational symmetry
        int : the number of duplicates

    Examples:
    >>> periodic_duplicates_in_cell([0,0,0])
    (([0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [1.0, 1.0, 0], [0, 0, 1.0], [1.0, 0, 1.0], [0, 1.0, 1.0], [1.0, 1.0, 1.0]), 8)
    >>> periodic_duplicates_in_cell([0,0.4,0])
    (([0, 0.4, 0], [1.0, 0.4, 0], [0, 0.4, 1.0], [1.0, 0.4, 1.0]), 4)
    '''
    from copy import deepcopy
    _pos = np.array(directCoord, dtype="float64")
    assert np.shape(_pos) == (3,)
    assert all(_pos - 1.0 < 0)
    _dupcs = []
    _dupcs.append(directCoord)
    # non-zero component
    _n = 2 ** (3 - np.count_nonzero(_pos))
    _iszero = _pos == 0
    for i in range(3):
        if _iszero[i]:
            _trans = deepcopy(_dupcs)
            for _c in _trans:
                _c[i] = 1.0
            _dupcs.extend(_trans)
    return tuple(_dupcs), _n
