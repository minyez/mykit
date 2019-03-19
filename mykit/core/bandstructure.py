# -*- coding: utf-8 -*-
'''
Module that defines class and utilities for band structure
'''

import numpy as np

from mykit.core.log import Verbose
from mykit.core.numeric import Prec

# eigen, occ (array): shape (nspins, nkpt, nbands)
DIM_EIGEN_OCC = 3
# pWave (array): shape (nspins, nkpt, nbands, natoms, nprojs)
KEYS_KPOINTS = ("coordinates", "weights")
KEYS_PROJECT = ("atoms", "projs", "pWave")
DIM_PROJECT = 5


class BandStructureError(Exception):
    pass


class BandStructure(Prec, Verbose):
    '''Base class for analyzing band structure data.

    The eigenvalues and occupations should be parsed in a shape of 
    (nspins, nkpts, nbands).
    The dimensions are automatically checked.
    Exception will be raised if their shapes are inconsistent.

    Optional keyword arguments, if parsed, should have consistent
    shape with the required arguments.
    Exception will be raised if their shapes are inconsistent.
    
    Args:
        eigen (array-like) : the energy (eigenvalues) of all bands to be considered
        occ (array-like) : the occupation numbers of all bands

    Optional args:
        kmesh (dict)
        projected (dict) : has three keys, "atoms, "projs" and "pWave"

    TODO:
        check kmesh consistency
    '''
    def __init__(self, eigen, occ, kmesh=None, projected=None):
        try:
            self._nspins, self._nkpts, self._nbands = \
                _check_eigen_occ_consistency(eigen, occ)
        except ValueError:
            info = "Bad eigen and occ shapes: {} vs {}".format(*map(np.shape, (eigen, occ)))
            raise BandStructureError(info)
        self._eigen = np.array(eigen, dtype=self._dtype)
        self._occ = np.array(occ, dtype=self._dtype)

        if projected != None:
            try:
                self._atoms, self._projs, pWave = \
                    self.__check_project_consistency(projected)
                self._pWave = np.array(pWave, dtype=self._dtype)
            except ValueError:
                self.print_warn("Bad projection input. Skip.")
        
    @property
    def eigen(self):
        return self._eigen
    @property
    def occ(self):
        return self._occ
    @property
    def atoms(self):
        if hasattr(self, "_atoms"):
            return self._atoms
        return None
    @property
    def projs(self):
        if hasattr(self, "_projs"):
            return self._projs
        return None
    @property
    def pWave(self):
        if hasattr(self, "_pWave"):
            return self._pWave
        return None
    @property
    def hasProjection(self):
        if self.pWave != None:
            return True
        return False
    
    def __check_project_consistency(self, projected):
        empty = ()
        try:
            assert isinstance(projected, dict)
        except AssertionError:
            return empty
        try:
            for key in KEYS_PROJECT:
                assert key in projected
        except AssertionError:
            return empty
        # print(projected)
        atoms = projected["atoms"]
        projs = projected["projs"]
        pWave = projected["pWave"]
        try:
            natoms = len(atoms)
            nprojs = len(projs)
            self.print_warn(np.shape(pWave), \
                (self._nspins, self._nkpts, self._nbands, natoms, nprojs),\
                level=3)
            assert np.shape(pWave) == \
                (self._nspins, self._nkpts, self._nbands, natoms, nprojs)
        except (AssertionError, TypeError):
            return empty
        return atoms, projs, pWave

    # def __check_kmesh_consistency(self, kmesh):
    #     try:
    #         for key in KEYS_KPOINTS:
    #             assert key in kpoints
    #     except AssertionError:
    #         return None


def _check_eigen_occ_consistency(eigen, occ):
    '''Check if eigenvalues and occupation number data have the same shape

    Both should have a shape 

    Returns:
        tuple, i.e. the shapes of eigen and occ are the same
    '''
    shapeEigen = np.shape(eigen)
    shapeOcc = np.shape(occ)
    if len(shapeEigen) == DIM_EIGEN_OCC and \
         len(shapeOcc) == DIM_EIGEN_OCC and \
            shapeEigen == shapeOcc:
        return shapeEigen
    return ()
