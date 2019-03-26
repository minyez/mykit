# -*- coding: utf-8 -*-
'''
Module that defines class and utilities for band structure
'''
import re
from numbers import Real

import numpy as np

from mykit.core.kmesh import check_kvecs_form_kpath
from mykit.core.log import Verbose
from mykit.core.numeric import Prec
from mykit.core.unit import EnergyUnit
from mykit.core.utils import get_str_indices_by_iden

# from mykit.core.visualizer import _BandVisualizer

# eigen, occ (array): shape (nspins, nkpt, nbands)
DIM_EIGEN_OCC = 3
'''int. Required dimension of eigenvalues and occupation numbers
'''

KEYS_BAND_PROJ = ("atoms", "projs", "pWave")
'''tuple. Required keys for wave projection input
'''

BAND_STR_PATTERN = re.compile(r"[vc]bm([-+][\d]+)?")

# DIM_PWAVE = 5
# '''pWave (array): shape (nspins, nkpt, nbands, natoms, nprojs)'''


class BandStructureError(Exception):
    pass


class BandStructure(Prec, Verbose, EnergyUnit):
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
        weight (array-like) : the weights of each kpoint, 1-d array,
            either int or float.

    Optional args:
        unit ('ev','ry','au'): the unit of the eigenvalues, in lower case. 
        kvec (array-like): kpoints vectors in reciprocal space (N.B. not the coordinate)
        efermi (float): the Fermi level. 
            If not parsed, the valence band maximum will be used.
        projected (dict) : wave projection information. 
            It should have three keys, "atoms", "projs" and "pWave"

    Attributes:

    '''

    def __init__(self, eigen, occ, weight, unit='ev', kvec=None, efermi=None, projected=None):
        try:
            self._nspins, self._nkpts, self._nbands = \
                _check_eigen_occ_weight_consistency(eigen, occ, weight)
        except ValueError:
            info = "Bad eigen, occ and weight shapes: {}, {}, {}".format(
                *map(np.shape, (eigen, occ, weight)))
            raise BandStructureError(info)
        try:
            self._eigen = np.array(eigen, dtype=self._dtype)
            self._occ = np.array(occ, dtype=self._dtype)
            self._weight = np.array(weight, dtype=self._dtype)
        except TypeError as _err:
            raise BandStructureError(_err)

        EnergyUnit.__init__(self, eunit=unit)
        self._emulti = {1: 2, 2: 1}[self._nspins]
        # channel: ispin, ikpt
        self._nelectPerChannel = np.sum(self._occ, axis=2) * self._emulti
        # In manual mode, one may parse the kpoints with all zero weight for band calculation
        # Reassign with unit weight
        if np.isclose(np.sum(self._weight), 0.0):
            self._weight[:] = 1.0
        self._nelectPerSpin = np.dot(
            self._nelectPerChannel, self._weight) / np.sum(self._weight)
        self._nelect = np.sum(self._nelectPerSpin)

        self._hasInftyCbm = False
        self._compute_vbm_cbm()
        if efermi is not None:
            assert isinstance(efermi, Real)
            self._efermi = efermi
        else:
            self._efermi = self.vbm

        self._projParsed = False
        self._parse_proj(projected)
        self._kvecParsed = False
        self._parse_kvec(kvec)

    def get_band_indices(self, *bands):
        '''Filter the band indices in ``bands``. 

        If no indices is specified, return all available band indices.

        Args:
            bands (int or str): the band identifier.
                can be band indices (int), or strings like "cbm", "vbm",
                "cbm-5", "vbm+2", etc.
        
        Returns:
            list
        '''
        if len(bands) == 0:
            b = list(range(self._nbands))
        else:
            b = []
            for ib in bands:
                if isinstance(ib, int):
                    if abs(ib) < self._nbands:
                        b.append(ib)
                elif isinstance(ib, str):
                    i = self._convert_band_str(ib)
                    if i is not None:
                        b.append(i)
        return b

    def _convert_band_str(self, bandStr):
        '''convert a string of band identifier, like "vbm", "cbm-2", "vbm+3"
        to the correpsonding band index.

        Args:
            bandStr (str)
        
        Returns:
            int
        '''
        assert isinstance(bandStr, str)
        if re.match(BAND_STR_PATTERN, bandStr):
            ref = {"v": self.ivbm[-1], "c": self.icbm[-1]}[bandStr[0]]
            if len(bandStr) == 3:
                return ref
            n = int(re.split(r"[-+]", bandStr)[-1])
            if bandStr[3] == '-':
                ib = ref - n
            else:
                ib = ref + n
            # check if the band index is valid
            if 0 <= ib < self.nbands:
                return ib
        return None

    @property
    def unit(self):
        return self._eunit

    @unit.setter
    def unit(self, newu):
        coef = self._get_eunit_conversion(newu)
        toConv = [self._eigen, self._bandWidth, 
                  self._vbmPerSpin, self._vbmPerChannel,
                  self._cbmPerSpin, self._cbmPerChannel,
                  ]
        if coef != 1:
            self._efermi *= coef
            self._vbm *= coef
            self._cbm *= coef
            for item in toConv:
                item *= coef
            self._eunit = newu.lower()

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
    def efermi(self):
        '''float. The energy of Fermi level.'''
        return self._efermi

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

        None if no projection information was parsed
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
        self._isKpath = None
        if not kvec is None:
            try:
                assert np.shape(kvec) == (self._nkpts, 3)
                self._kvec = np.array(kvec, dtype=self._dtype)
                # clean up very small values
                ind = self._kvec < 1.0E-10  # where values are low
                self._kvec[ind] = 0.0
                self._kvecParsed = True
                self._isKpath = False
                kSegments = check_kvecs_form_kpath(self._kvec)
                if kSegments != []:
                    self._isKpath = True
                    self._kLineSegs = kSegments
                else:
                    # use the first and last kpoints for plotting
                    self._kLineSegs = [(0, self.nkpts - 1),]
            except (AssertionError, ValueError):
                self.print_warn("Bad kpoint vectors input. Skip")

    @property
    def hasKvec(self):
        '''Bool'''
        return self._kvecParsed

    @property
    def isKpath(self):
        '''Bool. If the parsed kpoint vectors form a path in reciprocal space.

        None if kvec was not parsed.
        '''
        return self._isKpath

    @property
    def kvec(self):
        '''Array. kpoints vectors.

        None if kvec was not parsed.
        '''
        if self.hasKvec:
            return self._kvec
        return None
    
    @property
    def kLineSegs(self):
        '''List. Each member a list, containing the indices of starting and end
        kpoint vector
        '''
        return self._kLineSegs 

    def _generate_kpath_x(self):
        '''Generate the x coordinates for plotting a kpath

        x for the ending point of one line segment is the same to
        that for the starting point of next line segment.

        Returns
            list, each member a list of floats
        '''
        xs = []
        if not (self._isKpath and self.hasKvec):
            xs.append(list(range(self._nkpts)))
        else:
            _accuL = 0.0
            for sti, edi in self._kLineSegs:
                l = np.linalg.norm(self._kvec[sti, :] - self._kvec[edi, :])
                x = [_accuL + ik * l/(edi-sti) for ik in range(edi-sti+1)]
                xs.append(x)
                _accuL += l
        return xs

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
        # when any two indices of ivbm differ, the system is metal
        if np.max(self._ivbmPerChannel) == np.min(self._ivbmPerChannel):
            self._isMetal = False
        else:
            self._isMetal = True
        self._icbmPerChannel = self._ivbmPerChannel + 1
        # avoid IndexError when ivbm is the last band by imposing icbm = ivbm in this case
        ivbIsLast = self._ivbmPerChannel == self.nbands - 1
        if np.any(ivbIsLast):
            self.print_warn(
                "nbands {} is too small to get CB".format(self.nbands))
            self._icbmPerChannel[ivbIsLast] = self.nbands - 1

        self._vbmPerChannel = np.zeros(
            (self.nspins, self.nkpts), dtype=self._dtype)
        self._cbmPerChannel = np.zeros(
            (self.nspins, self.nkpts), dtype=self._dtype)
        # ? maybe optimize
        for i in range(self.nspins):
            for j in range(self.nkpts):
                vb = self._ivbmPerChannel[i, j]
                self._vbmPerChannel[i, j] = self.eigen[i, j, vb]
                if vb == self.nbands - 1:
                    self._hasInftyCbm = True
                    info = "VBM index for spin-kpt channel ({},{}) equals nbands.".format(
                        i+1, j+1)
                    self.print_warn(
                        info, "CBM for this channel set to infinity")
                    self._cbmPerChannel[i, j] = np.infty
                else:
                    self._cbmPerChannel[i, j] = self.eigen[i, j, vb+1]
        self._bandWidth = np.zeros(
            (self.nspins, self.nbands, 2), dtype=self._dtype)
        self._bandWidth[:, :, 0] = np.min(self._eigen, axis=1)
        self._bandWidth[:, :, 1] = np.max(self._eigen, axis=1)
        # VB indices
        self._ivbmPerSpin = np.array(((0, 0),)*self.nspins)
        self._vbmPerSpin = np.max(self._vbmPerChannel, axis=1)
        self._ivbmPerSpin[:, 0] = np.argmax(self._vbmPerChannel, axis=1)
        for i in range(self.nspins):
            ik = int(self._ivbmPerSpin[i, 0])
            self._ivbmPerSpin[i, 1] = self._ivbmPerChannel[i, ik]
        self._ivbm = np.array((0, 0, 0))
        self._ivbm[0] = int(np.argmax(self._vbmPerSpin))
        self._ivbm[1:3] = self._ivbmPerSpin[self._ivbm[0], :]
        self._vbm = self._vbmPerSpin[self._ivbm[0]]
        # CB indices
        self._icbmPerSpin = np.array(((0, 0),)*self.nspins)
        self._cbmPerSpin = np.min(self._cbmPerChannel, axis=1)
        self._icbmPerSpin[:, 0] = np.argmin(self._cbmPerChannel, axis=1)
        for i in range(self.nspins):
            ik = int(self._icbmPerSpin[i, 0])
            self._icbmPerSpin[i, 1] = self._icbmPerChannel[i, ik]
        self._icbm = np.array((0, 0, 0))
        self._icbm[0] = int(np.argmin(self._cbmPerSpin))
        self._icbm[1:3] = self._icbmPerSpin[self._icbm[0], :]
        self._cbm = self._cbmPerSpin[self._icbm[0]]

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

        int, shape (nspins, 2), ikpt, iband
        '''
        return self._ivbmPerSpin

    @property
    def icbmPerSpin(self):
        '''indices of conduction band minimum per spin

        int, shape (nspins, 2), ikpt, iband
        '''
        return self._icbmPerSpin

    @property
    def ivbm(self):
        '''index of valence band maximum

        int, shape (3,), ispin, ikpt, iband
        '''
        return self._ivbm

    @property
    def icbm(self):
        '''index of conduction band minimum

        int, shape (3,), ispin, ikpt, iband
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
    def bandWidth(self):
        '''the lower and upper bound of a band

        float, shape (nspins, nbands, 2)
        '''
        return self._bandWidth

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
            for key in KEYS_BAND_PROJ:
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
            self.print_log("Shapes: ", np.shape(pWave),
                           (self._nspins, self._nkpts,
                            self._nbands, natoms, nprojs),
                           level=3)
            assert np.shape(pWave) == \
                (self._nspins, self._nkpts, self._nbands, natoms, nprojs)
        except (AssertionError, TypeError):
            return ()
        return atoms, projs, pWave

    # * Projection related functions
    def effective_gap(self, ivb=None, atomVbm=None, projVbm=None,
                      icb=None, atomCbm=None, projCbm=None):
        '''Compute the effective band gap between ``ivb`` and ``icb``, 
        the responsible transition of which associates projector `projVbm` on `atomVbm` in VB
        and `projCbm` on atom `atomCbm` in CB.

        If no projection information was parsed, the inverse of the k-averaged gap inverse
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
        vbCoeffs = self.sum_atom_proj_comp(atomVbm, projVbm, fail_one=True)
        cbCoeffs = self.sum_atom_proj_comp(atomCbm, projCbm, fail_one=True)
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

    def sum_atom_proj_comp(self, atom, proj, fail_one=True):
        '''Sum the partial wave for projectors `proj` on atoms `atom`

        If the instance does not have projection information, 1 will be returned.

        Args:
            atom (int, str, Iterable)
            proj (int, str, Iterable)
            fail_one (bool): control the return when no projection is available. 
                if set True, return np.ones with correct shape, otherwise np.zeros

        Returns:
            (nspins, nkpts, nbands)
        '''
        if not self.hasProjection:
            func = {True: np.ones, False: np.zeros}
            try:
                return func[fail_one]((self.nspins, self.nkpts, self.nbands))
            except KeyError:
                raise TypeError("fail_one should be bool type.")
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
