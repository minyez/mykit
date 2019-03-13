# -*- coding: utf-8 -*-
'''Module that defines classes for crystal cell manipulation and symmtetry operation

The ``cell`` class and its subclasses accept the following kwargs when being instantialized:

    - unit (str): The unit system  to use, either "ang" (default) or "au".
    - coordSys (str): Coordinate system for the internal positions, 
      either "D" (Direct, default) or "C" (Cartesian)
    - allRelax (bool) : default selective dynamics option for atoms. 
      Set True (default) to allow all DOFs to relax
    - selectDyn (dict) : a dictionary with key-value pair as ``int: [bool, bool, bool]``, which controls
      the selective dynamic option for atom with the particular index (starting from 0). Default is an
      empty ``dict``
    - comment (str): message about the cell
'''
from collections import OrderedDict
# import spglib
from numbers import Real

import numpy as np

from mykit.core.constants import ang2au, au2ang, pi
from mykit.core.log import verbose
from mykit.core.numeric import prec


# ==================== classes ====================
class CellError(Exception):
    '''Exception in cell module
    '''
    pass


class Cell(prec, verbose):
    '''Cell structure class

    Args:
        latt (array-like) : The lattice vectors
        atoms (list of str) : The list of strings of type for each atom corresponding to the member in pos
        pos (array-like) : The internal coordinates of atoms

    Note:
        Various kwargs are acceptable for ``Cell`` and its subclasses. Check ``cell`` module docstring.
    
    Examples:
    >>> Cell([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], ["C"], [[0.0, 0.0, 0.0]])
    <mykit.core.cell.Cell>
    '''
    def __init__(self, latt, atoms, pos, **kwargs):

        # if kwargs:
            # super(lattice, self).__init__(**kwargs)
        self.comment= 'Default Cell class'
        self.__allRelax = True
        self.__selectDyn = {}
        self.__unit = 'ang'
        self.__coordSys = 'D'

        try:
            self.__latt = np.array(latt, dtype=self._dtype)
            self.__pos = np.array(pos, dtype=self._dtype)
        except ValueError as _err:
            raise CellError("Fail to create cell and pos array. Please check.")
        self.__atoms = [_a.capitalize() for _a in atoms]
        self.__parse_cellkw(**kwargs)
        # check input consistency
        self.__check_consistency()
        # move all atoms into the lattice (0,0,0)
        self.__assure_atoms_in_fist_lattice()
        # sanitize atoms arrangement
        self.__sanitize_atoms()

    def __len__(self):
        return len(self.__atoms)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.__pos[index, :]
        if isinstance(index, str):
            return self.get_sym_index(index)

    def __str__(self):
        return "{}\nLattice:\n{}\nAtoms: {}\nPositions:\n{}\nUnit: {}\nCoordinates in {}".format(\
            self.comment, self.latt, self.atoms, self.pos, self.unit, self.coordSys)
    
    def __repr__(self):
        return '{}'.format({"comment": self.comment, "latt": self.__latt, "atoms": self.__atoms, \
            "pos": self.__pos, "unit": self.__unit, "coordSys": self.__coordSys, \
            "allRelax": self.__allRelax, "selectDyn": self.__selectDyn})

    def __parse_cellkw(self, **kwargs):
        # accept_kw = ['unit', 'coordSys', 'allRelax', 'selectDyn']
        # for kw, v in kwargs.items():
        #     if kw in accept_kw:
        #         if kw == 'unit':
        #             self.__setattr__("__"+kw, v.lower())
        #             continue
        #         if kw == 'coordSys':
        #             self.__setattr__("__"+kw, v.upper())
        #             continue
        #         self.__setattr__("__"+kw, v)
        if 'unit' in kwargs:
            self.__unit = kwargs['unit'].lower()
        if 'coordSys' in kwargs:
            self.__coordSys = kwargs['coordSys'].upper()
        if "allRelax" in kwargs:
            self.__allRelax = kwargs["allRelax"]
        if "selectDyn" in kwargs:
            self.__selectDyn = kwargs["selectDyn"]
        if "comment" in kwargs:
            self.comment = kwargs["comment"]

    def get_cell(self):
        '''Purge out the cell, atoms and pos.
        
        ``cell``, ``atoms`` and ``pos`` are the minimal for constructing ``Cell`` and its subclasses.
        They can also be used to build sysmetry operations with spglib utilities.
        '''
        return self.__latt, self.__atoms, self.__pos

    def get_kwargs(self):
        '''return all kwargs useful to constract program-dependent cell input from ``Cell`` instance

        Returns :
            dictionary that can be parsed to ``create_from_cell`` class method.
        '''
        _d = {
                "unit" : self.__unit, 
                "coordSys" : self.__coordSys,
                "comment": self.comment,
                "allRelax" : self.__allRelax,
                "selectDyn" : self.__selectDyn
             }
        return _d

    def __check_consistency(self):
        try:
            assert self.__coordSys in ["C", "D"]
            assert self.__unit in ["ang", "au"]
            assert np.shape(self.__latt) == (3, 3)
            assert self.natoms > 0
            assert np.shape(self.__pos) == (self.natoms, 3)
        except AssertionError:
            raise CellError("Invalid cell setup")
        # ? switch automatically, or let the user deal with it
        try:
            assert self.vol > 0
        except AssertionError:
            raise CellError("Left-handed system found (vol<0). Switch two vector.")

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
            raise CellError("Fail to switch two atoms with indices {} and {}".format(iat1, iat2))

        self.__pos[[iat1, iat2]] = self.__pos[[iat2, iat1]]
        self.__atoms[iat1], self.__atoms[iat2] = self.__atoms[iat2], self.__atoms[iat1]

        _sfd1 = self.__selectDyn.pop(iat1,[])
        _sfd2 = self.__selectDyn.pop(iat2,[])
        if _sfd1 != []:
            self.__selectDyn.update({iat2: _sfd1})
        if _sfd2 != []:
            self.__selectDyn.update({iat1: _sfd2})

    def __assure_atoms_in_fist_lattice(self):
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

        Note that this is equivalent to ``latt[csymbol]``, given ``Cell`` instance latt.

        Args:
            csymbol (str) : chemical-symbol-like identifier
        '''
        __ind = []
        _csym = csymbol.capitalize()
        for _i, _s in enumerate(self.atoms):
            if _csym == _s:
                __ind.append(_i)
        return __ind

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
        for _i in range(_n -1):
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
        If ``reverse`` is set as False, atom with higher coordinate in lattice (0,0,0)
        will appear earlier, otherwise later

        Args :
            axis (1,2,3)
            reverse (bool)
        '''
        assert axis in range(1,4)
        assert isinstance(reverse, bool)
        __sortKeys = self.pos[:, axis-1]
        # self.print_log("__sortKeys in sort_pos", __sortKeys, level=3)
        for _at in self.atomTypes:
            __ind = self.get_sym_index(_at)
            self.__bubble_sort_atoms(__sortKeys, __ind, reverse=not reverse)

    # * Cell manipulation
    def scale(self, scale):
        '''Scale the lattice, i.e. increase the lattice by ``scale`` time
        '''
        assert isinstance(scale, Real)
        self.__latt = self.__latt * scale
        if self.__coordSys == "C":
            self.__pos = self.__pos * scale

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
        _shift = np.array([0.5,0.5,0.5], dtype=self._dtype) - _center

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
        return self.__latt

    @property
    def alen(self):
        return np.array([np.linalg.norm(x) for x in self.__latt], dtype=self._dtype)

    @property
    def lattConsts(self):
        '''Lattice constant of the cell, i.e., a, b, c, alpha, beta, gamma (In degree)
        '''
        _alen = self.alen
        _angle = []
        for i in range(3):
            j = (i+1)%3
            k = (i+2)%3
            _cos = np.dot(self.__latt[j], self.__latt[k])/_alen[j]/_alen[k]
            _angle.append(np.arccos(_cos))
        # convert to degree
        _angle = np.array(_angle, dtype=self._dtype) / pi * 180.0
        return (*_alen, *_angle)

    @property
    def latt(self):
        return self.__latt

    @property
    def atoms(self):
        return self.__atoms

    @property
    def pos(self):
        return self.__pos

    @property
    def unit(self):
        return self.__unit
    @unit.setter
    def unit(self, u):
        _u = u.lower()
        if _u == self.__unit:
            return
        _convDict = {"ang": au2ang, "au": ang2au}
        _conv = _convDict.get(_u)
        if _conv is not None:
            if self.__coordSys == "C":
                self.__pos = self.__pos * _conv
            self.__latt = self.__latt * _conv
            self.__unit = _u
        else:
            raise CellError("the length unit can only be either 'ang' (Angstrom) or 'au' (Bohr).")

    @property
    def coordSys(self):
        return self.__coordSys
    @coordSys.setter
    def coordSys(self, s):
        _s = s.upper()
        _convDict = {"C": self.__latt, "D": np.linalg.inv(self.__latt)}
        _conv = _convDict.get(_s)
        if _conv is not None:
            self.__pos = np.matmul(self.__pos, _conv)
        else:
            raise CellError("Only support \"D\" direct or fractional and \"C\" Cartisian coordinate.")

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
        return np.linalg.det(self.__latt)
    
    @property
    def natoms(self):
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
            b.append(np.cross(self.latt[j,:], self.latt[k,:]))
        return np.array(b, dtype=self._dtype) / self.vol

    @property
    def bIn2Pi(self):
        return self.recpLattIn2Pi

    @property
    def recpLatt(self):
        '''Reciprocal lattice vectors in unit^-1
        '''
        return self.recpLattIn2Pi * 2.0E0 * pi
    
    @property
    def b(self):
        '''Alias of reciprocal lattice vector
        '''
        return self.recpLatt

    @property
    def blen(self):
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
        _new = {}
        if len(iats) == 0:
            pass
        else:
            for _ia in iats:
                if _ia in range(self.natoms):
                    _new.update({_ia: select_dyn_flag_from_axis(axis, relax=False)})
            self.__set_sdFlags(_new)

    def set_relax(self, *iats, axis=0):
        '''Relax the atoms with index in iats

        Args:
            iats (list of int): the indices of atoms to relax
            axis (int or list): the axes along which the position of atom is relaxed
                It can be 0|1|2|3, or a list with all its members 1|2|3
        '''
        _new = {}
        if len(iats) == 0:
            pass
        else:
            for _ia in iats:
                if _ia in range(self.natoms):
                    _new.update({_ia: select_dyn_flag_from_axis(axis, relax=True)})
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
        assert isinstance(selectDyn, dict)
        for _k in selectDyn:
            _flag = selectDyn[_k]
            try:
                assert isinstance(_flag, list)
                assert len(_flag) == 3
                assert all([isinstance(_x, bool) for _x in _flag])
            except AssertionError:
                raise CellError("Bad flag for selective dynamics")
            else:
                self.__selectDyn.update({_k: _flag})

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
            _flag = [self.__allRelax,] *3
        else:
            _flag = [[self.__allRelax,]*3 for _r in range(self.natoms)]
            for _i in self.__selectDyn:
                _flag[_i] = self.__selectDyn[_i]
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
        '''Initialize a ``cell`` instance from a JSON file
        '''
        import json
        import os
        if pathJson is None or not os.path.isfile(pathJson):
            raise CellError("JSON file not found: {}".format(pathJson))
        with open(pathJson, 'r') as h:
            try:
                _j = json.load(h)
            except json.JSONDecodeError:
                raise CellError("invalid JSON file for cell: {}".format(pathJson))
        
        _argList = ["latt", "atoms", "pos"]
        _args = []
        for _i, arg in enumerate(_argList):
            v = _j.pop(arg, None)
            if v is None:
                raise CellError("invalid JSON file for cell: {}. No {}".format(pathJson, arg))
            _args.append(v)
        return cls(*_args, **_j)

    @classmethod
    def _bravais_o(cls, kind, atom, a, b, c, **kwargs):
        assert kind in ["P", "I", "F"]
        _latt = [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
        if kind == "P":
            _atoms = [atom,]
            _pos = [[0.0, 0.0, 0.0]]
        if kind == "I":
            _atoms = [atom,]*2
            _pos = [[0.0, 0.0, 0.0],[0.5,0.5,0.5]]
        if kind == "F":
            _atoms = [atom,]*4
            _pos = [[0.0, 0.0, 0.0],[0.0,0.5,0.5],[0.5,0.0,0.5],[0.5,0.5,0.0]]
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
        if not "comment" in kwargs:
            kwargs.update({"comment": "Simple orthorhombic lattice {}".format(atom)})
        return cls._bravais_o("P", atom, a, b, c, **kwargs)

    @classmethod
    def bravais_oI(cls, atom, a=1.0, b=2.0, c=3.0, **kwargs):
        '''Generate a body-centered orthorhombic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            a,b,c (float) : the lattice constants (a,b,c)
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        if not "comment" in kwargs:
            kwargs.update({"comment": "Body-centered orthorhombic lattice {}".format(atom)})
        return cls._bravais_o("I", atom, a, b, c, **kwargs)

    @classmethod
    def bravais_oF(cls, atom, a=1.0, b=2.0, c=3.0, **kwargs):
        '''Generate a face-centered orthorhombic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            a,b,c (float) : the lattice constants (a,b,c)
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        if not "comment" in kwargs:
            kwargs.update({"comment": "Face-centered orthorhombic lattice {}".format(atom)})
        return cls._bravais_o("F", atom, a, b, c, **kwargs)

    @classmethod
    def bravais_cP(cls, atom, aLatt=1.0, **kwargs):
        '''Generate a simple cubic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(aLatt)
        _latt = [[_a,0.0,0.0],[0.0,_a,0.0],[0.0,0.0,_a]]
        _atoms =[atom]
        _pos = [[0.0,0.0,0.0]]
        kwargs.pop("coordSys", None)
        if not "comment" in kwargs:
            kwargs.update({"comment": "Simple cubic lattice {}".format(atom)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def bravais_cI(cls, atom, aLatt=1.0, primitive=False, **kwargs):
        '''Generate a body-centered cubic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            primitive (bool) : if set True, the primitive cell will be generated
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(aLatt)
        if primitive:
            _latt = [[-_a/2.0,_a/2.0,_a/2.0],[_a/2.0,-_a/2.0,_a/2.0],[_a/2.0,_a/2.0,-_a/2.0]]
            _atoms =[atom]
            _pos = [[0.0,0.0,0.0]]
        else:
            _latt = [[_a,0.0,0.0],[0.0,_a,0.0],[0.0,0.0,_a]]
            _atoms =[atom,]*2
            _pos = [[0.0,0.0,0.0],[0.5,0.5,0.5]]
        kwargs.pop("coordSys", None)
        if not "comment" in kwargs:
            kwargs.update({"comment": "BCC {}".format(atom)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def bravais_cF(cls, atom, aLatt=1.0, primitive=False, **kwargs):
        '''Generate a face-centered cubic Bravais lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            primitive (bool) : if set True, the primitive cell will be generated
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(aLatt)
        if primitive:
            _latt = [[0.0,_a/2.0,_a/2.0],[_a/2.0,0.0,_a/2.0],[_a/2.0,_a/2.0,0.0]]
            _atoms =[atom]
            _pos = [[0.0,0.0,0.0]]
        else:
            _latt = [[_a,0.0,0.0],[0.0,_a,0.0],[0.0,0.0,_a]]
            _atoms =[atom,]*4
            _pos = [[0.0,0.0,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5],[0.5,0.5,0.0]]
        kwargs.pop("coordSys", None)
        if not "comment" in kwargs:
            kwargs.update({"comment": "FCC {}".format(atom)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def zincblende(cls, atom1, atom2, aLatt=1.0, primitive=False, **kwargs):
        '''Generate a standardized zincblende lattice (space group 216)

        ``atom1`` are placed at vertex and ``atom2`` at tetrahedron interstitial

        Args:
            atom1 (str): symbol of atom at vertex
            atom2 (str): symbol of atom at tetrahedron interstitial
            aLatt (float): the lattice constant of the conventional cell.
            primitive (bool): if set True, the primitive cell will be generated.
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(aLatt)
        if primitive:
            _latt = [[0.0, _a/2.0, _a/2.0], [_a/2.0, 0.0, _a/2.0], [_a/2.0, _a/2.0, 0.0]]
            _atoms = [atom1, atom2]
            _pos = [[0.0, 0.0, 0.0],
                    [0.25, 0.25, 0.25]]
        else:
            _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _a]]
            _atoms = [atom1,]*4 + [atom2,]*4
            _pos = [[0.0, 0.0, 0.0],
                    [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5],
                    [0.5, 0.5, 0.0],
                    [0.25, 0.25, 0.25],
                    [0.25, 0.75, 0.75],
                    [0.75, 0.25, 0.75],
                    [0.75, 0.75, 0.25]]
        kwargs.pop("coordSys", None)
        if not "comment" in kwargs:
            kwargs.update({"comment": "Zincblende {}{}".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)

    @classmethod
    def diamond(cls, atom, aLatt=1.0, primitive=False, **kwargs):
        '''Generate a standardized diamond lattice (space group 227)

        Args:
            atom (str): symbol of atom
            aLatt (float): the lattice constant of the conventional cell.
            primitive (bool): if set True, the primitive cell will be generated.
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        return cls.zincblende(atom, atom, aLatt=aLatt, primitive=primitive, **kwargs)

    @classmethod
    def wurtzite(cls, atom1, atom2, a=1.0, **kwargs):
        '''Generate a standardized wurtzite lattice (space group 186)

        ``atom1`` are placed at vertex and ``atom2`` at tetrahedron interstitial

        Args:
            atom1 (str): symbol of atom at vertex of lattice
            atom2 (str): symbol of atom at edge of lattice
            a (float): the lattice constant of the cell.
            kwargs: keyword argument for ``Cell`` except ``coordSys``
        '''
        _a = abs(a)
        _halfa = _a/2.0
        _c = _a * np.sqrt(8.0/3)
        _latt = [[_a, 0.0, 0.0], [-_halfa, np.sqrt(3)*_halfa, 0.0], [0.0, 0.0, _c]]
        _atoms = [atom1,]*2 + [atom2,]*2
        _pos = [[0.0, 0.0, 0.0],
                [2.0/3, 1.0/3, 0.5],
                [0.0, 0.0, 2.0/3],
                [2.0/3, 1.0/3, 1.0/6]]
        kwargs.pop("coordSys", None)
        if not "comment" in kwargs:
            kwargs.update({"comment": "Wurtzite {}{}".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)
    
    @classmethod
    def rutile(cls, atom1, atom2, a=1.0, c=2.0, u=0.31, **kwargs):
        '''Generate a standardized rutile lattice (space group 136)

        ``atom1`` are placed at vertex and ``atom2`` at tetrahedron interstitial

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
            raise CellError("Internal coordinate should be in (0,1), get {}".format(u))
        _a = abs(a)
        _c = abs(c)
        _latt = [[_a, 0.0, 0.0], [0.0, _a, 0.0], [0.0, 0.0, _c]]
        _atoms = [atom1,]*2 + [atom2,]*4
        _pos = [[0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [ u, u, 0.0],
                [-u,-u, 0.0],
                [0.5-u,0.5+u,0.5],
                [0.5+u,0.5-u,0.5]]
        kwargs.pop("coordSys", None)
        if not "comment" in kwargs:
            kwargs.update({"comment": "Rutile {}{}2".format(atom1, atom2)})
        return cls(_latt, _atoms, _pos, **kwargs)

def atoms_from_sym_nat(syms, nats):
    '''Generate ``atom`` list for ``Cell`` initilization from list of atomic symbols and number of atoms

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
        _list.extend([_s,] * _n)
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
        if axis in range(1,4):
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
                if _a in range(1,4):
                    _aList.append(_a)
    return tuple(_aList)


def periodic_duplicates_in_cell(directCoord):
    '''Return the coordinates and numbers of the duplicates of an atom in a cell due to lattice translation symmetry

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
