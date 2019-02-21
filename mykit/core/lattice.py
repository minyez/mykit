# -*- coding: utf-8 -*-
'''Module that defines classes for lattice manipulation and symmtetry operation
'''
import numpy as np
# import spglib
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
    
    Available ``kwargs``:
        unit (str): The system of unit to use. Either "ang" or "au"
        coordSys (str): Coordinate system for the internal positions. Either "D" (Direct) or "C" (Cartesian)
        fSelectDyn (bool) : switch of selective dynamics for geometry optimization
        selectDyn (dict) : 
    
    Examples:
    >>> latt = lattice([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], ["C"], [[0.0, 0.0, 0.0]])
    '''
    def __init__(self, cell, atoms, pos, **kwargs):

        # if kwargs:
            # super(lattice, self).__init__(**kwargs)
        self.comment= ''
        self._fSelectDyn = False
        self._selectDyn = None
        self._natoms = 0
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

    def __parse_kwargs(self, **kwargs):
        if 'unit' in kwargs:
            self.__unit = kwargs['unit'].lower()
        if 'coordSys' in kwargs:
            self.__coordSys = kwargs['coordSys'].upper()
        
        _directAssign = {
                            "comment": self.comment,
                            "fSelectDyn": self._fSelectDyn,
                            "selectDyn": self._selectDyn
                        }

        for _k in _directAssign:
            if _k in kwargs:
                _directAssign[_k] = kwargs[_k]


    def get_kwargs(self):
        '''return all kwargs useful to constract program-dependent lattice input from ``lattice`` instance

        Returns :
            dictionary that can be parsed to ``create_from_lattice`` class method.
        '''
        _d = {
                "unit" : self.__unit, 
                "coordSys" : self.__coordSys,
                "comment": self.comment,
                "fSelectDyn" : self._fSelectDyn,
                "selectDyn" : self._selectDyn
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

    def move(self, iAtom):
        pass

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

    def __len__(self):
        return len(self.__atoms)

    def __getitem__(self, index):
        return self.__pos[index, :]

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

    # def rebuild(self):
    #     '''Rebuild the lattice from cell, atoms and pos
    #     '''
    #     pass


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
    
    Examples:
    >>> generate_atoms_from_sym_nat(["C", "Al", "F"], [2, 3, 1])
    ["C", "C", "Al", "Al", "Al", "F"]
    '''
    assert len(syms) == len(nats)
    _list = []
    for _s, _n in zip(syms, nats):
        _list.extend([_s,] * _n)
    return _list