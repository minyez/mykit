# -*- coding: utf-8 -*-
'''Module that defines classes for lattice manipulation and symmtetry operation
'''
import numpy as np
# import spglib
from mykit.core.log import verbose
from mykit.core.numeric import prec
from mykit.core.constants import pi, au2ang, ang2au

_latt_kwargs_doc = '''
'''

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

    Available ``kwargs``:
        unit (str): The system of unit to use. Either "ang" or "au"
        coordSys (str): Coordinate system for the internal positions. Either "D" (Direct) or "C" (Cartesian)
        allRelax (bool) : default selective dynamics option for atoms. Set True to allow all DOFs to relax
        selectDyn (dict) : a dictionary with key-value pair as "int: [bool, bool, bool]", which controls
            the selective dynamic option for atom with the particular index (starting from 0)
    
    Examples:
    >>> lattice([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], ["C"], [[0.0, 0.0, 0.0]])
    <mykit.core.lattice.lattice>
    '''
    def __init__(self, cell, atoms, pos, **kwargs):

        # if kwargs:
            # super(lattice, self).__init__(**kwargs)
        self.comment= ''
        self.__allRelax = True
        self.__selectDyn = {}
        self.__unit = 'ang'
        self.__coordSys = 'D'

        try:
            self.__cell = np.array(cell, dtype=self._dtype)
            self.__pos = np.array(pos, dtype=self._dtype)
        except ValueError as _err:
            raise latticeError("Fail to create cell and pos array. Please check.")
        self.__atoms = atoms
        self.__parse_kwargs(**kwargs)
        # check input consistency
        self.__check_consistency()
        # self.__generate_index_dict()

    def __len__(self):
        return len(self.__atoms)

    def __getitem__(self, index):
        return self.__pos[index, :]

    def __parse_kwargs(self, **kwargs):
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

    def _switch_two_atom_index(self, iat1, iat2):
        '''switch the index of atoms with index iat1 and iat2

        Except ``__pos`` and ``__atoms``, 
        this method should also deals possible switch in other positional 
        attributes, e.g.

        Note that this method is mainly for sorting, and does NOT change 
        the geometry of the lattice at all.
        
        - ``selectDyn`` (DONE)
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

    def __sanitize(self):
        '''Sanitize the lattice properties after creation. 
        '''
        raise NotImplementedError    
    
    def __sort(self):
        '''Sort the atoms by its coordinate
        '''
        raise NotImplementedError    

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

    def get_latt(self):
        '''Purge the cell, atoms and pos, which is minimal for constructing ``lattice`` and its subclasses
        '''
        return self.__cell, self.__atoms, self.__pos

    def __check_consistency(self):
        try:
            assert self.__coordSys in ["C", "D"]
            assert self.__unit in ["ang", "au"]
            assert np.shape(self.__cell) == (3, 3)
            assert np.shape(self.__pos) == (len(self.__atoms), 3)
        except AssertionError:
            raise latticeError("Invalid lattice setup")

    # TODO move atom
    def move(self, ia):
        raise NotImplementedError

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
        return tuple(_list)

    @property
    def typeIndex(self):
        _ats = self.atomTypes
        _dict = {}
        for i, _at in enumerate(_ats):
            _dict.update({_at: i})
        return tuple([ _dict[_a] for _a in self.__atoms])

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

    def set_fix(self, *iats, axis=0):
        '''Fix the atoms with index in iats

        Args:
            axis (int or list): the axis along which the position of atom is fixed
                It can be 0|1|2|3, or a list with all its members 1|2|3
        '''
        _new = {}
        if len(iats) == 0:
            pass
        else:
            for _ia in iats:
                if _ia in range(self.natoms):
                    _new.update({_ia: select_dyn_flag_from_axis(axis, relax=True)})
            self.set_sdFlags(_new)

    def set_relax(self, *iats, axis=0):
        '''Relax the atoms with index in iats
        '''
        _new = {}
        if len(iats) == 0:
            pass
        else:
            for _ia in iats:
                if _ia in range(self.natoms):
                    _new.update({_ia: select_dyn_flag_from_axis(axis, relax=False)})
            self.set_sdFlags(_new)

    def set_sdFlags(self, selectDyn):
        assert isinstance(selectDyn, dict)
        for _k in selectDyn:
            _flag = selectDyn[_k]
            try:
                assert isinstance(_flag, list)
                assert len(_flag) == 3
                assert all([isinstance(_x, bool) for _x in _flag])
            except AssertionError:
                pass
            else:
                self.__selectDyn.update({_k: _flag})

    def sdFlags(self, ia=-1):
        '''Return the selective dynamic flag (bool)

        Args:
            ia (int) : index of atom

        Returns:
            3-member list, if ia is in range(natoms),
            otherwise natoms-member list, each member a 3-member list
            as the flag for that atom
        '''
        self.print_log("Global selective dynamics flag: {}".format(self.__allRelax), level=3)
        if ia in self.__selectDyn:
            self.print_log("Found custom flag in __selectDyn for atom {}".format(ia), depth=1, level=3)
            _flag = self.__selectDyn[ia]
        elif ia in range(self.natoms):
            self.print_log("Use global flag for atom {}".format(ia), level=3, depth=1)
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
#         vol = np.linalg.det(self.a)
#         self.vol = vol

#         # reciprocal lattice vectors in 2pi
#         b = []
#         for i in range(3):
#             j = (i + 1) % 3
#             k = (i + 2) % 3
#             b.append(np.cross(self.a[j], self.a[k]))
#         b = np.array(b, dtype=self._dtype)
#         self.b2Pi = np.divide(b, vol)
#         self.b = np.multiply(b, 2.0 * pi)

#         # length of lattice vectors
#         self.bLen = np.array([np.linalg.norm(x) for x in self.b], dtype=self._dtype)
#         self.bLen2Pi = np.array([np.linalg.norm(x) for x in self.b2Pi], dtype=self._dtype)
#         self.volBZ2Pi = np.linalg.det(self.b2Pi)
#         self.volBZ = np.linalg.det(b)

    # Factory methods
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
        return cls.__bravis_c(atom, 'P', aLatt=aLatt, **kwargs)

    @classmethod
    def bravis_cI(cls, atom, aLatt=1.0, **kwargs):
        '''Generate a body-centered cubic Bravis lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            unit (str) : the unit system to use
        '''
        return cls.__bravis_c(atom, 'I', aLatt=aLatt, **kwargs)

    @classmethod
    def bravis_cF(cls, atom, aLatt=1.0, **kwargs):
        '''Generate a face-centered cubic Bravis lattice

        Args:
            atom (str) : the chemical symbol of atom
            alatt (float) : the lattice constant (a)
            unit (str) : the unit system to use
        '''
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
    >>> generate_atoms_from_sym_nat(["C", "Al", "F"], [2, 3, 1])
    ["C", "C", "Al", "Al", "Al", "F"]
    '''
    assert len(syms) == len(nats)
    _list = []
    for _s, _n in zip(syms, nats):
        _list.extend([_s,] * _n)
    return _list


def sym_nat_from_atoms(atoms):
    '''Generate lists of atomic symbols and number of atoms

    The order of appearence of the element is conserved in the output.

    Args :
        atoms (list of str) : symbols of each atom in the lattice

    Returns :
        list of str : atomic symbols
        list of int : number of atoms for each symbol
    
    Examples:
    >>> generate_atoms_from_sym_nat(["C", "Al", "Al", "C", "Al", "F"])
    ["C", "Al", "F"], [2, 3, 1]
    '''
    _syms = []
    _natsDict = {}
    for _at in atoms:
        if _at in _syms:
            _natsDict[_at] += 1
        else:
            _syms.append(_at)
            _natsDict.update({_at: 0})
    return _syms, [_natsDict[_at] for _at in _syms]


def select_dyn_flag_from_axis(axis, relax=True):
    assert isinstance(relax, bool)
    _base = [relax, relax, relax]
    if isinstance(axis, int):
        if axis == 0:
            _base = [not _b for _b in _base]
        if axis in range(1,4):
            _base[axis-1] = not _base[axis-1]
    elif isinstance(axis, (list, tuple)):
        __axisSet = list(set(axis))
        for _a in __axisSet:
            assert isinstance(_a, int)
            if _a == 0:
                _base = [not _b for _b in _base]
                break
            if _a in range(1,4):
                _base[_a-1] = not _base[_a-1]
    return _base