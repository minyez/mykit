# -*- coding: utf-8 -*-
'''Module that defines classes for lattice manipulation and symmtetry operation

The ``lattice`` class and its subclasses accept the following kwargs when being instantialized:

    - unit (str): The unit system  to use, either "ang" (default) or "au".
    - coordSys (str): Coordinate system for the internal positions, 
      either "D" (Direct, default) or "C" (Cartesian)
    - allRelax (bool) : default selective dynamics option for atoms. 
      Set True (default) to allow all DOFs to relax
    - selectDyn (dict) : a dictionary with key-value pair as ``int: [bool, bool, bool]``, which controls
      the selective dynamic option for atom with the particular index (starting from 0). Default is an
      empty ``dict``
    - comment (str): message about the lattice
'''
import numpy as np
# import spglib
from numbers import Real
from mykit.core.log import verbose
from mykit.core.numeric import prec
from mykit.core.constants import pi, au2ang, ang2au

# ==================== classes ====================
class latticeError(Exception):
    '''Exception in lattice module
    '''
    pass


class lattice(prec, verbose):
    '''Lattice structure class

    Args:
        cell (array-like) : The lattice vectors
        atoms (list of str) : The list of strings of type for each atom corresponding to the member in pos
        pos (array-like) : The internal coordinates of atoms

    Note:
        Various kwargs are acceptable for ``lattice`` and its subclasses. Check ``lattice`` module docstring.
    
    Examples:
    >>> lattice([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], ["C"], [[0.0, 0.0, 0.0]])
    <mykit.core.lattice.lattice>
    '''
    def __init__(self, cell, atoms, pos, **kwargs):

        # if kwargs:
            # super(lattice, self).__init__(**kwargs)
        self.comment= 'Default lattice class'
        self.__allRelax = True
        self.__selectDyn = {}
        self.__unit = 'ang'
        self.__coordSys = 'D'

        try:
            self.__cell = np.array(cell, dtype=self._dtype)
            self.__pos = np.array(pos, dtype=self._dtype)
        except ValueError as _err:
            raise latticeError("Fail to create cell and pos array. Please check.")
        self.__atoms = [_a.capitalize() for _a in atoms]
        self.__parse_lattkw(**kwargs)
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
        return "{}\nCell:\n {}\nAtoms: {}\nPositions:\n {}\nUnit: {}\nCoordinates in {}".format(self.comment, self.cell, self.atoms, self.pos, self.unit, self.coordSys)
    
    def __repr__(self):
        return self.__str__()

    def __parse_lattkw(self, **kwargs):
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

    def get_latt(self):
        '''Purge out the cell, atoms and pos.
        
        ``cell``, ``atoms`` and ``pos`` are the minimal for constructing ``lattice`` and its subclasses.
        They can also be used to build sysmetry operations with spglib utilities.
        '''
        return self.__cell, self.__atoms, self.__pos

    def get_kwargs(self):
        '''return all kwargs useful to constract program-dependent lattice input from ``lattice`` instance

        Returns :
            dictionary that can be parsed to ``create_from_lattice`` class method.
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
            assert np.shape(self.__cell) == (3, 3)
            assert self.natoms > 0
            assert np.shape(self.__pos) == (self.natoms, 3)
        except AssertionError:
            raise latticeError("Invalid lattice setup")
        # ? switch automatically, or let the user deal with it
        try:
            assert self.vol > 0
        except AssertionError:
            raise latticeError("Left-handed system found (vol<0). Switch two vector.")

    def _switch_two_atom_index(self, iat1, iat2):
        '''switch the index of atoms with index iat1 and iat2

        Except ``__pos`` and ``__atoms``, 
        this method should also deals possible switch in other positional 
        attributes, e.g.
        
        - ``selectDyn`` (DONE)

        Note that this method is mainly for sorting use, and does NOT change 
        the geometry of the lattice at all.
        '''
        try:
            assert iat1 in range(self.natoms)
            assert iat2 in range(self.natoms)
            assert iat1 != iat2
        except AssertionError:
            raise latticeError("Fail to switch two atoms with indices {} and {}".format(iat1, iat2))

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
        '''Get the indices of atom with element symbol ``csymbol``

        Note that this is equivalent to ``latt[csymbol]``, given ``lattice`` instance latt.

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
        '''Scale the lattice, i.e. increase the lattice cell by ``scale`` time
        '''
        assert isinstance(scale, Real)
        self.__cell = self.__cell * scale
        if self.__coordSys == "C":
            self.__pos = self.__pos * scale

    # TODO add more atoms
    def add_atom(self):
        raise NotImplementedError

    # TODO move atom
    def __move(self, ia):
        raise NotImplementedError

    @property
    def center(self):
        '''Calculate the center of all atoms in the cell
        '''
        assert self.coordSys == "D"
        _posSum = np.zeros(3, dtype=self._dtype)
        _n = 0
        for i in range(self.natoms):
            _dupcs, _dn = periodic_duplicates_in_cell(self.__pos[i, :])
            _n += _dn
            for _dupc in _dupcs:
                np.add(_posSum, _dupc, _posSum)
        # _posSum = np.sum(self.__pos, axis=0)
        # check periodic duplicate, by recognizing number of zeros in pos.
        # _dup = 3 - np.count_nonzero(self.__pos, axis=1)
        # _dup = np.power(2, _dup)
        # _n = np.sum(_dup)
        # _posSum = np.sum(self.__pos * _dup[:, None], axis=0)
        return _posSum / _n

    # ! Not work now
    def centering(self, axis=0):
        '''Centering the atoms along axes. Mainly use for slab model.

        Args:
            axis (int or iterable of int) : the axes along which the atoms will be centered.
        '''
        _aList = axis_list(axis)
        _wasDirect = self.coordSys == "D"
        if not _wasDirect:
            self.coordSys = "D"
        # get the geometric center of all atoms
        _center = self.center
        _shift = np.array([0.5,0.5,0.5], dtype=self._dtype) - _center

        for i in range(3):
            ia = i + 1
            if ia not in _aList:
                _shift[i] = 0.0
        self.__pos = np.add(self.__pos, _shift)
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
        if not _wasDirect:
            self.coordSys = "C"

    @property
    def a(self):
        return self.__cell

    @property
    def alen(self):
        return np.array([np.linalg.norm(x) for x in self.__cell], dtype=self._dtype)

    @property
    def lattConsts(self):
        '''Lattice constant of the lattice, i.e., a, b, c, alpha, beta, gamma (In degree)
        '''
        _alen = self.alen
        _angle = []
        for i in range(3):
            j = (i+1)%3
            k = (i+2)%3
            _cos = np.dot(self.__cell[j], self.__cell[k])/_alen[j]/_alen[k]
            _angle.append(np.arccos(_cos))
        # convert to degree
        _angle = np.array(_angle, dtype=self._dtype) / pi * 180.0
        return (*_alen, *_angle)

    @property
    def cell(self):
        return self.__cell

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
        _convDict = {"ang": au2ang, "au": ang2au}
        _conv = _convDict.get(_u)
        if _conv is not None:
            if self.__coordSys == "C":
                self.__pos = self.__pos * _conv
            self.__cell = self.__cell * _conv
        else:
            raise latticeError("the length unit can only be either 'ang' (Angstrom) or 'au' (Bohr).")

    @property
    def coordSys(self):
        return self.__coordSys
    @coordSys.setter
    def coordSys(self, s):
        _s = s.upper()
        _convDict = {"C": self.__cell, "D": np.linalg.inv(self.__cell)}
        _conv = _convDict.get(_s)
        if _conv is not None:
            self.__pos = np.matmul(self.__pos, _conv)

    @property
    def atomTypes(self):
        _list = []
        for _a in self.__atoms:
            if _a not in _list:
                _list.append(_a)
        return _list

    @property
    def typeIndex(self):
        _ats = self.atomTypes
        _dict = {}
        for i, _at in enumerate(_ats):
            _dict.update({_at: i})
        return [_dict[_a] for _a in self.__atoms]

    @property
    def vol(self):
        return np.linalg.det(self.__cell)
    
    @property
    def natoms(self):
        return len(self.__atoms)

    @property
    def useSelDyn(self):
        if self.__allRelax and not bool(self.__selectDyn):
            return False
        return True

    @property
    def recpCellIn2Pi(self):
        '''Reciprocal lattice vectors in 2Pi unit^-1
        '''
        b = []
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            b.append(np.cross(self.cell[j,:], self.cell[k,:]))
        return np.array(b, dtype=self._dtype) / self.vol

    @property
    def bIn2Pi(self):
        return self.recpCellIn2Pi

    @property
    def recpCell(self):
        '''Reciprocal lattice vectors in unit^-1
        '''
        return self.recpCellIn2Pi * 2.0E0 * pi
    
    @property
    def b(self):
        '''Alias of reciprocal lattice vector
        '''
        return self.recpCell

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

    def __set_sdFlags(self, selectDyn):
        assert isinstance(selectDyn, dict)
        for _k in selectDyn:
            _flag = selectDyn[_k]
            try:
                assert isinstance(_flag, list)
                assert len(_flag) == 3
                assert all([isinstance(_x, bool) for _x in _flag])
            except AssertionError:
                raise latticeError("Bad flag for selective dynamics")
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

# TODO make reciprocal lattice information to properties
#    def __calc_latt(self):
#         '''Calculate lattice information from the input
#         '''

#         b = np.array(b, dtype=self._dtype)
#         self.b2Pi = np.divide(b, vol)
#         self.b = np.multiply(b, 2.0 * pi)
#         # length of lattice vectors
#         _b = np.array([np.linalg.norm(x) for x in self.b], dtype=self._dtype)
#         self.volBZ2Pi = np.linalg.det(self.b2Pi)
#         self.volBZ = np.linalg.det(b)

    # * Factory methods
    @classmethod
    def __bravis_c(cls, atom, typeLatt, aLatt=1.0, **kwargs):

        _type = typeLatt.upper()
        try:
            assert isinstance(aLatt, (int, float))
        except AssertionError:
            raise latticeError("alatt should be a number")
        try:
            assert _type in ["P", "I", "F"]
        except AssertionError:
            raise latticeError("Invalid cubic Bravis system")

        _a = abs(aLatt)
        _cell = [[_a,0.0,0.0],[0.0,_a,0.0],[0.0,0.0,_a]]
        if _type == "P":
            _atoms =[atom]
            _pos = [[0.0,0.0,0.0]]
        if _type == "I":
            _atoms =[atom, atom]
            _pos = [[0.0,0.0,0.0],[0.5,0.5,0.5]]
        if _type == "F":
            _atoms =[atom, atom, atom, atom]
            _pos = [[0.0,0.0,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5],[0.5,0.5,0.0]]

        return cls(_cell, _atoms, _pos, **kwargs)

    @classmethod
    def bravis_cP(cls, atom, aLatt=1.0, **kwargs):
        '''Generate a primitive cubic Bravis lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            unit (str) : the unit system to use
        '''
        if not "comment" in kwargs:
            kwargs.update({"comment": "Primitive cubic lattice"})
        return cls.__bravis_c(atom, 'P', aLatt=aLatt, **kwargs)

    @classmethod
    def bravis_cI(cls, atom, aLatt=1.0, **kwargs):
        '''Generate a body-centered cubic Bravis lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            unit (str) : the unit system to use
        '''
        if not "comment" in kwargs:
            kwargs.update({"comment": "Body-centered cubic lattice"})
        return cls.__bravis_c(atom, 'I', aLatt=aLatt, **kwargs)

    @classmethod
    def bravis_cF(cls, atom, aLatt=1.0, **kwargs):
        '''Generate a face-centered cubic Bravis lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            unit (str) : the unit system to use
        '''
        if not "comment" in kwargs:
            kwargs.update({"comment": "Face-centered cubic lattice"})
        return cls.__bravis_c(atom, 'F', aLatt=aLatt, **kwargs)
    

# class symmetry(prec):
#     '''Class to obtain symmetry of a lattice

#     Args:
#         latt (lattice): the :class:`lattice` instance for symmetry search
#     '''
#     def __init__(self, latt):
#         pass


def atoms_from_sym_nat(syms, nats):
    '''Generate ``atom`` list for ``lattice`` initilization from list of atomic symbols and number of atoms

    Args :
        syms (list of str) : atomic symbols
        nats (list of int) : number of atoms for each symbol

    Returns :
        a list of str, containing symbol of each atom in the lattice
    
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
        atoms (list of str) : symbols of each atom in the lattice

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

# TODO general this function to mirrors in n-th lattice shell
def periodic_duplicates_in_cell(directCoord):
    '''Return the coordinates and numbers of the duplicates of an atom in a cell due to lattice translation symmetry

    Args:
        directCoord (array): the direct coordinate of an atom in the cell
    
    Note:
        The function works only when each component belongs to [0,1)
    
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