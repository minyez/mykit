# -*- coding: utf-8 -*-
'''
Module that defines class and utilities for band structure
'''

from numbers import Real

import numpy as np

from mykit.core.log import Verbose
from mykit.core.numeric import Prec
from mykit.core.utils import get_str_indices_by_iden

# eigen, occ (array): shape (nspins, nkpt, nbands)
DIM_EIGEN_OCC = 3
'''int. Required dimension of eigenvalues and occupation numbers
'''

KEYS_PROJ = ("atoms", "projs", "pWave")
'''tuple. Required keys for projection input
'''

# DIM_PWAVE = 5
# '''pWave (array): shape (nspins, nkpt, nbands, natoms, nprojs)'''


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
        weight (array-like) : the integer weights of each kpoint, 1-d int array

    Optional args:
        kvec (array-like): kpoints vectors in reciprocal space (N.B. not the coordinate)
        efermi (float): the Fermi level. If not parsed, the valence band maximum will be used.
        projected (dict) : partial wave information, which has three keys, "atoms, "projs" and "pWave"

    Attributes:

    '''
    def __init__(self, eigen, occ, weight, kvec=None, efermi=None, projected=None):
        try:
            self._nspins, self._nkpts, self._nbands = \
                _check_eigen_occ_weight_consistency(eigen, occ, weight)
        except ValueError:
            info = "Bad eigen, occ and weight shapes: {}, {}, {}".format(*map(np.shape, (eigen, occ, weight)))
            raise BandStructureError(info)
        try:
            self._eigen = np.array(eigen, dtype=self._dtype)
            self._occ = np.array(occ, dtype=self._dtype)
            self._weight = np.array(weight, dtype=self._dtype)
        except TypeError as _err:
            raise BandStructureError(_err)
        self._emulti = {1: 2, 2: 1}[self._nspins]
        # channel: ispin, ikpt
        self._nelectPerChannel = np.sum(self._occ, axis=2) * self._emulti
        self._nelectPerSpin = np.dot(self._nelectPerChannel, self._weight) / sum(self._weight)
        self._nelect = np.sum(self._nelectPerSpin)

        self._hasInftyCbm = False
        self._compute_vbm_cbm()
        if not efermi is None:
            assert isinstance(efermi, Real)
            self._efermi = efermi
        else:
            self._efermi = np.max(self.vbm)
    
        self._projParsed = False
        self._parse_proj(projected)
        self._kvecParsed = False
        self._parse_kvec(kvec)
    
    @property
    def eigen(self):
        '''Array. Eigenvalues, (nspins, nkpts, nbands)'''
        return self._eigen
    @property
    def occ(self):
        '''Array. Occupation numbers, (nspins, nkpts, nbands)'''
        return self._occ
    @property
    def nelect(self):
        '''Float. Total number of electrons'''
        return self._nelect
    @property
    def nspins(self):
        '''Int. number of spin channels'''
        return self._nspins
    @property
    def nkpts(self):
        '''Int. number of kpoints'''
        return self._nkpts
    @property
    def nbands(self):
        '''Int. number of bands'''
        return self._nbands

    def _parse_proj(self, projected):
        '''Parse the partial wave information
        '''
        if not projected is None:
            try:
                self._atoms, self._projs, pWave = \
                    self._check_project_consistency(projected)
                self._pWave = np.array(pWave, dtype=self._dtype)
                self._projParsed = True
            except ValueError:
                self.print_warn("Bad projection input. Skip.")

    @property
    def hasProjection(self):
        '''Bool'''
        return self._projParsed
    @property
    def atoms(self):
        '''list of strings, types of each atom. 
        
        None if no projection information was parsed
        '''
        if self.hasProjection:
            return self._atoms
        return None
    @property
    def natoms(self):
        try:
            return len(self.atoms)
        except TypeError:
            return 0
    @property
    def projs(self):
        '''list of strings, names of atomic projectors
        
        None if no projection information wa parsed
        '''
        if self.hasProjection:
            return self._projs
        return None
    @property
    def nprojs(self):
        try:
            return len(self.projs)
        except TypeError:
            return 0
    @property
    def pWave(self):
        '''Array, partial waves for each projector on each atom.

        shape (nspins, nkpts, nbands, natoms, nprojs)

        None if no projection information was parsed.
        '''
        if self.hasProjection:
            return self._pWave
        return None

    def _parse_kvec(self, kvec):
        '''Parse the kpoints vectors
        '''
        if not kvec is None:
            try:
                assert np.shape(kvec) == (self._nkpts, 3)
                self._kvec = np.array(kvec, dtype=self._dtype)
                # clean up very small values
                ind = self._kvec < 1.0E-10  # where values are low
                self._kvec[ind] = 0.0
                self._kvecParsed = True
            except (AssertionError, ValueError):
                self.print_warn("Bad kpoint vectors input. Skip")
    @property
    def hasKvec(self):
        '''Bool'''
        return self._kvecParsed
    @property
    def kvec(self):
        '''Array. kpoints vectors.
        
        None if kvec was not parsed.
        '''
        if self.hasKvec:
            return self._kvec
        return None

    def _compute_vbm_cbm(self):
        '''Get the indices of band edges for each spin-kpoint channel

        Note:
            When nbands is too small and all bands are found to be valence bands,
            or it is the case for one spin-kpoint channel, the corresponding CB will
            be set to ``np.infty``. 
            Setup of indices of CB remains, thus IndexError might be raised when trying
            get CB value from `icbm` attributes
        '''
        isEmpty = self._occ < (self._thresEmp * self._emulti)
        self._ivbmPerChannel = np.sum(isEmpty, axis=2) - 1
        self._icbmPerChannel = self._ivbmPerChannel + 1
        self._vbmPerChannel = np.zeros((self.nspins, self.nkpts), dtype=self._dtype)
        self._cbmPerChannel = np.zeros((self.nspins, self.nkpts), dtype=self._dtype)

        # ? maybe optimize
        for i in range(self.nspins):
            for j in range(self.nkpts):
                vb = self._ivbmPerChannel[i, j]
                self._vbmPerChannel[i, j] = self.eigen[i, j, vb]
                if vb == self.nbands - 1:
                    self._hasInftyCbm = True
                    info = "VBM index for spin-kpt channel ({},{}) equals nbands.".format(i+1, j+1)
                    self.print_warn(info, "CBM for this channel set to infinity")
                    self._cbmPerChannel[i, j] = np.infty
                else:
                    self._cbmPerChannel[i, j] = self.eigen[i, j, vb+1]

        self._ivbmPerSpin = np.max(self._ivbmPerChannel, axis=1)
        self._vbmPerSpin = np.max(self._vbmPerChannel, axis=1)
        self._ivbm = np.max(self._ivbmPerSpin)
        self._vbm = np.max(self._vbmPerSpin)

        self._icbmPerSpin = self._ivbmPerSpin + 1
        self._cbmPerSpin = np.min(self._cbmPerChannel, axis=1)
        self._icbm = self._ivbm + 1
        self._cbm = np.min(self._cbmPerSpin)

    @property
    def ivbm(self):
        '''indices of valence band maximum at each spin-kpt channel
        
        int, shape (nspins, nkpts)
        '''
        return self._ivbmPerChannel
    @property
    def icbm(self):
        '''indices of conduction band minimum at each spin-kpt channel
        
        int, shape (nspins, nkpts)
        '''
        return self._icbmPerChannel
    @property
    def vbm(self):
        '''valiues of valence band maximum at each spin-kpt channel
        
        float, shape (nspins, nkpts)
        '''
        return self._vbmPerChannel
    @property
    def cbm(self):
        '''values of conduction band minimum at each spin-kpt channel
        
        float, shape (nspins, nkpts)
        '''
        return self._cbmPerChannel
    @property
    def directGap(self):
        '''Direct gap at each spin-kpt channel

        float, shape (nspins, nkpts)
        '''
        return self._cbmPerChannel - self._vbmPerChannel
    @property
    def fundGap(self):
        '''Fundamental gap for each spin channel

        float, shape (nspins,)
        '''
        return self._cbmPerSpin - self._vbmPerSpin
    @property
    def fundTrans(self):
        '''Transition responsible for the fundamental gap

        int, shape (nspins, 2)
        '''
        vb = np.argmax(self._vbmPerChannel, axis=1)
        cb = np.argmin(self._vbmPerChannel, axis=1)
        return tuple(zip(vb, cb))
    @property
    def kAvgGap(self):
        '''direct band gap averaged over kpoints

        float, shape (nspins,)
        '''
        return np.dot(self.directGap, self._weight) / np.sum(self._weight)

    def _check_project_consistency(self, projected):
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

    # * Projection related functions
    def effective_gap(self, atomVbm=None, projVbm=None, atomCbm=None, projCbm=None):
        '''Compute the effective band gap, 
        the responsible transition of which associates projector `projVbm` on `atomVbm` in VB
        and `projCbm` on atom `atomCbm` in CB

        See paper for more information

        Args:
            atomVbm (int, str, Iterable): atom where the VB projector is located
            atomCbm (int, str, Iterable): atom where the CB projector is located
            projVbm (int, str, Iterable): index of VB projector
            projCbm (int, str, Iterable): index of CB projector
        '''
        if not self.hasProjection:
            return None
        raise NotImplementedError

    def _sum_atom_proj_comp(self, atom, proj):
        '''Sum the partial wave for projector `proj` on atom `atom`

        Args:
            atom (int, str, Iterable):
            proj (int, str, Iterable):

        Returns:
            (nspins, nkpts, nbands)
        '''
        if not self.hasProjection:
            return 0.0
        

    def _get_atom_indices(self, atom):
        if not self.hasProjection:
            return []
        else:
            return get_str_indices_by_iden(self._atoms, atom)

    def _get_proj_indices(self, proj):
        if not self.hasProjection:
            return []
        else:
            return get_str_indices_by_iden(self._projs, proj)
        

def _check_eigen_occ_weight_consistency(eigen, occ, weight):
    '''Check if eigenvalues, occupation number and kweights data have the correct shape

    Returns:
        tuple, the shape of eigen/occ when the shapes of input are consistent,
            empty tuple otherwise.
    '''
    shapeEigen = np.shape(eigen)
    shapeOcc = np.shape(occ)
    shapeKw = np.shape(weight)
    consist = len(shapeEigen) == DIM_EIGEN_OCC and \
                len(shapeOcc) == DIM_EIGEN_OCC and \
                        shapeEigen == shapeOcc and \
                             len(shapeKw) == 1 and \
                  len(weight) == shapeEigen[1]
    if consist:
        return shapeEigen
    return ()
