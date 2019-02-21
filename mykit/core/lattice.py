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
        unit (str): The system of unit to use. Either "ang" or "au"
        coordSys (str): Coordinate system for the internal positions. Either "D" (Direct) or "C" (Cartesian)
    
    Examples:
    >>> latt = lattice([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]], ["C"], [[0.0, 0.0, 0.0]])
    '''
    def __init__(self, cell, atoms, pos, **kwargs):

        # if kwargs:
            # super(lattice, self).__init__(**kwargs)
        self.comment= ''
        self._atomIndexDict = {}
        self._natoms = 0
        self.__unit = 'ang'
        self.__coordSys = 'D'

        try:
            self.cell = np.array(cell, dtype=self._dtype)
            self.pos = np.array(pos, dtype=self._dtype)
        except ValueError as _err:
            raise latticeError("Fail to create cell and pos array. Please check.")
        self.atoms = atoms
        self.__parse_kwargs(**kwargs)

        # check input consistency
        self.__check_consistency()
        self.__generate_index_dict()

    def __parse_kwargs(self, **kwargs):
        if 'comment' in kwargs:
            self.comment = kwargs['comment']
        if 'unit' in kwargs:
            self.__unit = kwargs['unit'].lower()
        if 'coordSys' in kwargs:
            self.__coordSys = kwargs['coordSys'].upper()


    def get_kwargs(self):
        '''return all kwargs useful to constract program-dependent lattice input from ``lattice`` instance

        Returns :
            dictionary that can be parsed to ``create_from_lattice`` class method.
        '''
        _d = {
                "unit" : self.__unit, 
                "coordSys" : self.__coordSys,
                "comment": self.comment
             }
        return _d

    def __check_consistency(self):
        try:
            assert self.__coordSys in ["C", "D"]
            assert self.__unit in ["ang", "au"]
            assert np.shape(self.cell) == (3, 3)
            assert np.shape(self.pos) == (len(self.atoms), 3)
        except AssertionError:
            raise latticeError("Invalid lattice setup")

    def __generate_index_dict(self):
        '''Generate the atomType-index dictionary
        '''
        _aSet = tuple(set(self.atoms))
        for _i, _a in enumerate(_aSet):
            self._atomIndexDict.update({_a: _i})

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
                self.pos = self.pos * _conv
            self.cell = self.cell * _conv

    @property
    def coordSys(self):
        return self.__coordSys
    @coordSys.setter
    def coordSys(self, s):
        _s = s.upper()
        _convDict = {"C": self.cell, "D": np.linalg.inv(self.cell)}
        _conv = _convDict.get(_s)

        if _conv is not None:
            self.pos = np.matmul(self.pos, _conv)

    @property
    def vol(self):
        return np.linalg.det(self.cell) 
    
    @property
    def natoms(self):
        return len(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, index):
        return self.pos[index, :]

    # Factory methods
    @classmethod
    def __bravis_c(cls, atom, typeLatt, aLatt=1.0, **kwargs):

        __type = typeLatt.upper()
        try:
            assert isinstance(aLatt, (int, float))
        except AssertionError:
            raise latticeError("alatt should be a number")
        try:
            assert __type in ["P", "I", "F"]
        except AssertionError:
            raise latticeError("Invalid cubic Bravis system")

        _a = abs(aLatt)
        __cell = [[_a,0.0,0.0],[0.0,_a,0.0],[0.0,0.0,_a]]
        if __type == "P":
            __atoms =[atom]
            __pos = [[0.0,0.0,0.0]]
        if __type == "I":
            __atoms =[atom, atom]
            __pos = [[0.0,0.0,0.0],[0.5,0.5,0.5]]
        if __type == "F":
            __atoms =[atom, atom, atom, atom]
            __pos = [[0.0,0.0,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5],[0.5,0.5,0.0]]

        return cls(__cell, __atoms, __pos, **kwargs)

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