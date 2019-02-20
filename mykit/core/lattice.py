# -*- coding: utf-8 -*-
'''Module that defines classes for lattice manipulation and symmtetry operation
'''
import numpy as np
# import spglib
from mykit.core.numeric import prec
from mykit.core.constants import pi, au2ang

# ==================== classes ====================
class latticeError(Exception):
    '''Exception in lattice module
    '''
    pass


class lattice(prec):
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
    def __init__(self, cell, atoms, pos, unit="ang", coordSys="D", **kwargs):

        self._atomIndexDict = {}
        self._natoms = 0
        self._unit = ''
        self._coordSys = ''

        try:
            self.cell = np.array(cell, dtype=self._dtype)
            self.pos = np.array(pos, dtype=self._dtype)
        except ValueError as _err:
            raise latticeError("Fail to create cell and pos array. Please check.")
        self.atoms = atoms
        self._unit = unit.lower()
        self._coordSys = coordSys.upper()

        # check input consistency
        self.__check_consistency()
        self.__generate_index_dict()

    def __check_consistency(self):
        try:
            assert self._coordSys in ["C", "D"]
            assert self._unit in ["ang", "au"]
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
        return self._unit
    @unit.setter
    def unit(self, u):
        _u = u.lower()
        try:
            assert _u in ["ang", "au"]
        except AssertionError:
            pass
        else:
            if _u == self._unit:
                pass
            else:
                self._unit = _u
                if _u == "ang":
                    self.__au2ang()
                elif _u == "au":
                    self.__ang2au()

    def __au2ang(self):
        if self._coordSys == "C":
            self.pos = self.pos * au2ang
        self.cell = self.cell * au2ang

    def __ang2au(self):
        if self._coordSys == "C":
            self.pos = self.pos / au2ang
        self.cell = self.cell / au2ang

    @property
    def coordSys(self):
        return self._coordSys
    @coordSys.setter
    def coordSys(self, s):
        _s = s.upper()
        try:
            assert _s in ["D", "C"]
        except AssertionError:
            pass
        else:
            if _s == self._coordSys:
                pass
            else:
                self._coordSys = _s
                if _s == "D":
                    self.__cart2direct()
                elif _s == "C":
                    self.__direct2cart()
                    
    def __cart2direct(self):
        self.pos = np.matmul(self.pos, np.linalg.inv(self.cell))

    def __direct2cart(self):
        self.pos = np.matmul(self.pos, self.cell)


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
        
