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
    They will be ignored if their shapes are inconsistent.
    
    Args:
        eigen (array-like) : the energy (eigenvalues) of all bands to be considered
        occ (array-like) : the occupation numbers of all bands
        weight (array-like) : the integer weights of each kpoint, 1-d int array

    Optional args:
        kvec (array-like): kpoints vectors in reciprocal space (N.B. not the coordinate)
        efermi (float): the Fermi level. 
            If not parsed, the valence band maximum will be used.
        projected (dict) : partial wave information. 
            It should have three keys, "atoms, "projs" and "pWave"

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
        # In manual mode, one may parse the kpoints with all zero weight for band calculation
        # Reassign with unit weight
        if np.isclose(np.sum(self._weight), 0.0):
            self._weight[:] = 1.0
        self._nelectPerSpin = np.dot(self._nelectPerChannel, self._weight) / np.sum(self._weight)
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
    def weight(self):
        '''Array. kpoints weight (int like), (nkpts, )'''
        return self._weight
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
        '''compute the band edges on each spin-kpt channel

        Note:
            When nbands is too small and all bands are found to be valence bands,
            or it is the case for one spin-kpoint channel, the corresponding CB will
            be set to ``np.infty``. 
            Setup of indices of CB remains, thus IndexError might be raised when trying
            get CB value from `icbm` attributes
        '''
        isOcc = self._occ > self._thresOcc
        
        self._ivbmPerChannel = np.sum(isOcc, axis=2) - 1
        self._ivbmPerSpin = np.max(self._ivbmPerChannel, axis=1)
        self._ivbm = np.max(self._ivbmPerSpin)
        # when any two indices of ivbm differ, the system is metal
        if np.max(self._ivbmPerChannel) == np.min(self._ivbmPerChannel):
            self._isMetal = False
        else:
            self._isMetal = True

        # determine CB indices
        self._icbmPerChannel = self._ivbmPerChannel + 1
        # avoid IndexError when ivbm is the last band by imposing icbm = ivbm in this case
        ivbIsLast = self._ivbmPerChannel == self.nbands - 1
        if np.any(ivbIsLast):
            self.print_warn("nbands is too small to get CB")
            self._icbmPerChannel[ivbIsLast] = self.nbands - 1
        # self._icbmPerSpin = self._ivbmPerSpin + 1
        self._icbmPerSpin = np.min(self._icbmPerChannel, axis=1)
        # self._icbm = self._ivbm + 1
        self._icbm = np.min(self._icbmPerSpin)

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

        # self._vbmPerSpin = np.max(self._vbmPerChannel, axis=1)
        # self._cbmPerSpin = np.min(self._cbmPerChannel, axis=1)
        self._vbmPerSpin = np.zeros(self.nspins, dtype=self._dtype)
        self._cbmPerSpin = np.zeros(self.nspins, dtype=self._dtype)
        for i in range(self.nspins):
            self._vbmPerSpin[i] = np.max(self._eigen[i, :, self._ivbmPerSpin[i]])
            self._cbmPerSpin[i] = np.min(self._eigen[i, :, self._icbmPerSpin[i]])
        self._vbm = np.max(self._vbmPerSpin)
        self._cbm = np.min(self._cbmPerSpin)

    @property
    def isMetal(self):
        '''bool.
        
        True if the bandstructure belongs to a metal, False otherwise
        '''
        return self._isMetal
    @property
    def ivbmPerChannel(self):
        '''indices of valence band maximum at each spin-kpt channel
        
        int, shape (nspins, nkpts)
        '''
        return self._ivbmPerChannel
    @property
    def icbmPerChannel(self):
        '''indices of conduction band minimum at each spin-kpt channel
        
        int, shape (nspins, nkpts)
        '''
        return self._icbmPerChannel
    @property
    def ivbmPerSpin(self):
        '''indices of valence band maximum per spin
        
        int, shape (nspins,)
        '''
        return self._ivbmPerSpin
    @property
    def icbmPerSpin(self):
        '''indices of conduction band minimum per spin
        
        int, shape (nspins,)
        '''
        return self._icbmPerSpin
    @property
    def ivbm(self):
        '''index of valence band maximum
        
        int
        '''
        return self._ivbm
    @property
    def icbm(self):
        '''index of conduction band minimum
        
        int
        '''
        return self._icbm
    @property
    def vbmPerChannel(self):
        '''valiues of valence band maximum at each spin-kpt channel
        
        float, shape (nspins, nkpts)
        '''
        return self._vbmPerChannel
    @property
    def cbmPerChannel(self):
        '''values of conduction band minimum at each spin-kpt channel
        
        float, shape (nspins, nkpts)
        '''
        return self._cbmPerChannel
    @property
    def vbmPerSpin(self):
        '''valiues of valence band maximum per spin
        
        float, shape (nspins,)
        '''
        return self._vbmPerSpin
    @property
    def cbmPerSpin(self):
        '''values of conduction band minimum per spin
        
        float, shape (nspins,)
        '''
        return self._cbmPerSpin
    @property
    def vbm(self):
        '''value of valence band maximum
        
        float
        '''
        return self._vbm
    @property
    def cbm(self):
        '''value of conduction band minimum
        
        float
        '''
        return self._cbm
    @property
    def directGap(self):
        '''Direct gap between VBM and CBM of each spin-kpt channel

        float, shape (nspins, nkpts)
        '''
        return self._cbmPerChannel - self._vbmPerChannel
    @property
    def fundGap(self):
        '''Fundamental gap for each spin channel.
        
        float, shape (nspins,)
        If it is metal, it is equivalent to the negative value of bandwidth
        of the unfilled band.
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
    def effective_gap(self, ivb=None, atomVbm=None, projVbm=None, \
            icb=None, atomCbm=None, projCbm=None):
        '''Compute the effective band gap between ``ivb`` and ``icb``, 
        the responsible transition of which associates projector `projVbm` on `atomVbm` in VB
        and `projCbm` on atom `atomCbm` in CB.

        If projection information was not parsed, the inverse of the k-averaged gap inverse
        will be returned.

        Args:
            ivb (int): index of the lower band. Use VBM if not specified or is invalid index.
            icb (int): index of the upper band. Use CBM if not specified or is invalid index.
            atomVbm (int, str, Iterable): atom where the VB projector is located
            atomCbm (int, str, Iterable): atom where the CB projector is located
            projVbm (int, str, Iterable): index of VB projector
            projCbm (int, str, Iterable): index of CB projector

        Note:
            Spin-polarization is not considered in retriving projection coefficients.
        '''
        vbCoeffs = self._sum_atom_proj_comp(atomVbm, projVbm)
        cbCoeffs = self._sum_atom_proj_comp(atomCbm, projCbm)
        if ivb is None or not ivb in range(self.nbands):
            vbCoeff = vbCoeffs[:, :, np.max(self.ivbm)]
        else:
            vbCoeff = vbCoeffs[:, :, ivb]
        if icb is None or not icb in range(self.nbands):
            cbCoeff = cbCoeffs[:, :, np.min(self.icbm)]
        else:
            cbCoeff = cbCoeffs[:, :, icb]
        # ! abs is added in case ivb and icb are put in the opposite 
        inv = np.sum(np.abs(np.reciprocal(self.directGap) * vbCoeff * cbCoeff))
        if np.allclose(inv, 0.0):
            return np.infty
        return 1.0/inv

    def _sum_atom_proj_comp(self, atom, proj):
        '''Sum the partial wave for projectors `proj` on atoms `atom`

        If the instance does not have projection information, 1 will be returned.

        Args:
            atom (int, str, Iterable):
            proj (int, str, Iterable):

        Returns:
            (nspins, nkpts, nbands)
        '''
        if not self.hasProjection:
            return np.ones((self.nspins, self.nkpts, self.nbands))
        if atom is None:
            atInd = list(range(self.natoms))
        else:
            atInd = self._get_atom_indices(atom)
        if proj is None:
            prInd = list(range(self.nprojs))
        else:
            prInd = self._get_proj_indices(proj)
        coeff = np.zeros((self.nspins, self.nkpts, self.nbands))
        for a in atInd:
            for p in prInd:
                coeff += self.pWave[:, :, :, a, p]
        return coeff

    def _get_atom_indices(self, atom):
        if self.hasProjection:
            return get_str_indices_by_iden(self._atoms, atom)
        return []

    def _get_proj_indices(self, proj):
        if self.hasProjection:
            return get_str_indices_by_iden(self._projs, proj)
        return []
        

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
