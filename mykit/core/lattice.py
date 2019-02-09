# -*- coding: utf-8 -*-
'''Module that defines classes for lattice manipulation and symmtetry operation
'''
from math import pi
import numpy as np
import spglib
from mykit.core.numeric import prec
#try:
#    import ase
#except ImportError:
#    pass

# ==================== classes ====================
class lattice(prec):
    '''Lattice structure class

    Args:
        cell (array_like) : The lattice vectors
        atoms (list) : The list of strings of type for each atom corresponding to the member in pos
        pos (array_like) : The internal coordinates of atoms
        unit (str): The system of unit to use, either "ang" or "au"
        coordSys (str): either "D" (direct) or "C" (Cartesian)
    '''
    def __init__(self, cell, atoms, pos, unit="ang", coordSys="D"):
        assert coordSys.lower() in ["C", "D"]
        assert unit.lower() in ["ang", "au"]
        self.a = np.ndarray(cell, dtype=self.__dtype)
        self.atoms = atoms
        self.pos = np.ndarray(pos, dtype=self.__dtype)
        self.__atomIndexDict = {}
        self.atomIndex = []
        self.coordSys = coordSys
        self.unit = unit
        self.__generate_atomDict()
        self.__check_consistency()
        self.__calc_latt()

    def __generate_atom_index(self):
        __uAtomList = tuple(set(self.atoms))
        for i, __atom in enumerate(__uAtomList):
            self.__atomIndexDict.update({__atom: i})
        for __atom in self.atoms:
            self.atomIndex.append(self.__atomIndexDict[__atom])

    def __check_consistency(self):
        try:
            assert np.shape(self.a) == (3, 3)
            assert np.shape(self.pos) == (len(self.atoms), 3)
        except AssertionError:
            raise AssertionError("Inconsistent lattice import")

    def __calc_latt(self):
        '''
        '''
        vol = np.linalg.det(self.a)
        self.vol = vol

        # reciprocal lattice vectors in 2pi
        b = []
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            b.append(np.cross(self.a[j], self.a[k]))
        b = np.array(b, dtype=self.__dtype)
        self.b2Pi = np.divide(b, vol)
        self.b = np.multiply(b, 2.0 * pi)

        # length of lattice vectors
        self.bLen = np.array([np.linalg.norm(x) for x in self.b], dtype=self.__dtype)
        self.bLen2Pi = np.array([np.linalg.norm(x) for x in self.b2Pi], dtype=self.__dtype)
        self.volBZ2Pi = np.linalg.det(self.b2Pi)
        self.volBZ = np.linalg.det(b)

    def rebuild(self):
        '''Rebuild the lattice from cell, atoms and pos
        '''
        pass

class symmetry(prec):
    '''Class to obtain symmetry of a lattice

    Args:
        latt (lattice): the :class:`lattice` instance for symmetry search
    '''
    def __init__(self, latt):
        pass
        
