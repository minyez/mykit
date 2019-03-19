# -*- coding: utf-8 -*-
'''
Module that defines class and utilities for band structure
'''

from numbers import Real

import numpy as np

from mykit.core.log import Verbose
from mykit.core.numeric import Prec

# eigen, occ (array): shape (nspins, nkpt, nbands)
DIM_EIGEN_OCC = 3
# pWave (array): shape (nspins, nkpt, nbands, natoms, nprojs)
KEYS_KMESH = ("coordinates", "weights")
KEYS_PROJ = ("atoms", "projs", "pWave")
# DIM_PWAVE = 5


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
        efermi (float)
        kmesh (dict)
        projected (dict) : has three keys, "atoms, "projs" and "pWave"
        b (array-like): reciprocal lattice vectors with shape (3,3)

    TODO:
        check kmesh consistency
    '''
    def __init__(self, eigen, occ, efermi=None, kmesh=None, projected=None, b=None):
        try:
            self._nspins, self._nkpts, self._nbands = \
                _check_eigen_occ_consistency(eigen, occ)
        except ValueError:
            info = "Bad eigen and occ shapes: {} vs {}".format(*map(np.shape, (eigen, occ)))
            raise BandStructureError(info)
        self._eigen = np.array(eigen, dtype=self._dtype)
        self._occ = np.array(occ, dtype=self._dtype)
        self._emulti = {1: 2, 2: 1}[self._nspins]
        # channel: ispin, ikpt
        self._nelectPerChannel = np.sum(self._occ, axis=2) * self._emulti
        self._nelectPerSpin = np.sum(self._nelectPerChannel, axis=1) / self._nkpts
        self._nelect = np.sum(self._nelectPerSpin)

        self._vbm = None
        if not efermi is None:
            assert isinstance(efermi, Real)
            self._efermi = efermi
        else:
            self._efermi = self.vbm
    
        self._projParsed = False
        self.__parse_proj(projected)
        self._kmeshParsed = False
        self.__parse_kmesh(kmesh)
        self._bParsed = False
        self.__parse_b(b)
    
    def __parse_proj(self, projected):
        '''Parse the partial wave information
        '''
        if not projected is None:
            try:
                self._atoms, self._projs, pWave = \
                    self.__check_project_consistency(projected)
                self._pWave = np.array(pWave, dtype=self._dtype)
                self._projParsed = True
            except ValueError:
                self.print_warn("Bad projection input. Skip.")

    def __parse_kmesh(self, kmesh):
        '''Parse the kpoints mesh
        '''
        pass

    def __parse_b(self, b):
        '''Parse the reciprocal lattice vector
        '''
        if not b is None:
            try:
                assert np.shape(b) == (3,3)
                self._b = np.array(b, dtype=self._dtype)
                self._bParsed = True
            except (AssertionError, ValueError):
                self.print_warn("Bad reciprocal lattice vectors. Skip.")

    @property
    def eigen(self):
        return self._eigen
    @property
    def occ(self):
        return self._occ
    @property
    def nelect(self):
        return self._nelect
    @property
    def nspins(self):
        return self._nspins
    @property
    def nkpts(self):
        return self._nkpts
    @property
    def nbands(self):
        return self._nbands
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
        if not self.pWave is None:
            return True
        return False

    def __check_project_consistency(self, projected):
        try:
            assert isinstance(projected, dict)
        except AssertionError:
            return ()
        try:
            for key in KEYS_PROJ:
                assert key in projected
        except AssertionError:
            return ()
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
            return ()
        return atoms, projs, pWave

    @property
    def vbm(self):
        if self._vbm is None:
            return self.get_vbm()
        return self._vbm

    def get_vbm(self):
        vbm = 0.0
        return vbm

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
