# -*- coding: utf-8 -*-
'''
Module that defines class and utilities for band structure
'''
from abc import ABC, abstractmethod
from copy import deepcopy
from numbers import Real

import numpy as np

from mykit.core.kmesh import (KSYMBOL_LATEX, KmeshError,
                              check_kvecs_form_kpath, kpath_decoder)
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

# ==========================
# related with visualization
# ==========================

def init_band_visualizer(bs, plotter=None, **kwargs):
    '''The wrapper to return a Visualizer object for plotting a band structure.

    Args:
        bs (``BandStructure``): the band structure to visualize
        plotter (str): the plotter to use. Currently support 
            - "pyplot": matplotlib
        kwargs: keyword arguments to parse to initialize the plotter.
            - "pyplot": matplotlib.pyplot.subplots
    '''
    # Supported plotter. The first is default.
    _SUPPORT_PLOTTER = {
        'pyplot': BSVisualizerPyplot,
    }
    default = 'pyplot'
    try:
        assert isinstance(bs, BandStructure)
    except AssertionError:
        raise TypeError("BandStructure object should be parsed.")
    if plotter is None:
        p = default
    else:
        if plotter not in _SUPPORT_PLOTTER:
            raise NotImplementedError("plotter {} not supported.".format(plotter))
        else:
            p = plotter
    return _SUPPORT_PLOTTER[p](bs, **kwargs)


# contract for a visualizer of band structure
class _BSVisualizer(ABC, Verbose):
    '''Abstract base class for band structure visualizer
    with different plotter
    '''
    def __init__(self, bs, align_vbm=False):
        # check bandstructure
        assert isinstance(bs, BandStructure)
        if not bs.isKpath:
            self.print_warn("k-vectors of Band structure do not form a kpath. Plotting anyway...")
        self._bs = deepcopy(bs)
        self._xs = bs._generate_kpath_x()
        self._efermi = bs.efermi
        self._nspins = bs.nspins
        self._ispin = 0
        self.alignAtVbm = align_vbm
        self._drawnKsym = False
    
    @property
    def alignAtVbm(self):
        '''bool. whether the VBM is set as energy level'''
        return self._alignAtVbm
    
    @alignAtVbm.setter
    def alignAtVbm(self, newValue):
        if not isinstance(newValue, bool):
            raise TypeError("alignAtVbm should be bool")
        self._alignAtVbm = newValue

    @property
    def drawnSym(self):
        '''bool. If the kpoints symbols have been drawn'''
        return self._drawnKsym

    @property
    def ispin(self):
        '''int. index of spin channel to plot'''
        return self._ispin
    
    @ispin.setter
    def ispin(self, newValue):
        if not isinstance(newValue, int):
            raise TypeError("ispin should be int")
        elif newValue >= self._nspins:
            raise TypeError("spin channel index overflow. nspins = %d" % self._nspins)
        self._ispin = newValue
    
    # Abstract methods to implement
    @abstractmethod
    def set_title(self, title, **kwargs):
        pass

    @abstractmethod
    def set_elim(self, bottom, top, **kwargs):
        pass
        
    @abstractmethod
    def draw(self, *bands, **kwargs):
        pass

    @abstractmethod    
    def mark_ksymbols(self, kpath, **kwargs):
        pass
    
    @abstractmethod
    def draw_proj(self, atom, proj, *bands, **kwargs):
        pass
    @abstractmethod
    def export(self, *args, **kwargs):
        pass


class BSVisualizerPyplot(_BSVisualizer):
    '''The class to draw band structure by matplotlib.pyplot. subplots is used.

    Currently only support drawing one bandstructure, i.e. 1 Axes in Fig.
    Can only plot one spin channel.

    Args:
        bs (``BandStructure``): the band structure to plot
        kwargs: keyword arguments to parse to initialize with `pyplot.subplots`
    '''

    def __init__(self, bs, align_vbm=False, **kwargs):
        '''
        '''
        super(BSVisualizerPyplot, self).__init__(bs, align_vbm=align_vbm)
        try:
            import matplotlib.pyplot as plt
        except ModuleNotFoundError:
            raise ValueError("Matplotlib is required to use pyplot plotter")
        # initialize the Figure and Axes
        self._fig, self._axes = plt.subplots(**kwargs)
        # set xlimit
        self._axes.set_xlim([self._xs[0][0], self._xs[-1][-1]])
        # set x tick length to zero to make them invisible
        self._axes.tick_params(axis=u'x', which=u'both', length=0)
        self._axes.set_ylabel("Energy ({})".format(bs.unit))
    
    def set_title(self, title, **kwargs):
        '''Set the title of figure
        
        Wrapper to pyplot.Axes.set_title
        '''
        self._axes.set_title(title, **kwargs)

    def set_elim(self, bottom, top, **kwargs):
        '''Set the energy limit to show
        
        wrapper to Axes.set_ylim
        '''
        self._axes.set_ylim(bottom=bottom, top=top, **kwargs)

    def draw(self, *bands, **kwargs):
        '''draw the selected bands

        Args:
            bands (int): the indices of bands to plot
                All bands will be plot by default.
            kwargs: keywords to parse to Axes.plot
        '''
        # iterate for each line segments
        if 'color' not in kwargs:
            kwargs['color'] = 'k'
        bs = self._bs
        xs = self._xs
        if len(bands) == 0:
            b = list(range(bs.nbands))
        else:
            b = bands
        for i, (stk, edk) in enumerate(bs.kLineSegs):
            for ib in b:
                eigen = bs.eigen[self.ispin, stk:edk+1, ib] - \
                    bs.vbmPerSpin[self.ispin] * int(self.alignAtVbm)
                self._axes.plot(xs[i], eigen, **kwargs)
            # draw vertical line to separate the different line segments
            if i != len(bs.kLineSegs) - 1:
                self._axes.axvline(xs[i][-1], color="k")
        # draw Fermi level
        if self.alignAtVbm:
            self._axes.axhline(0.0, color="k", linestyle='dashed')
        else:
            self._axes.axhline(bs.efermi, color="k", linestyle='dashed')
        
    def mark_ksymbols(self, kpath):
        '''Mark the kpoints symbols on the plot

        Args
            kpath (str): the kpath string.
                If the kpoints of BandStructure do not form a kpath, skip.
                The string will be decoded by `kpath_decoder` to a list of kpoints symbols.
                followed by a consistency check.
                If the check fail, a warning will be prompted and no symbols plotted
        '''
        if not self._bs.isKpath:
            return
        try:
            ksyms = kpath_decoder(kpath)
        except KmeshError:
            return
        if len(ksyms) / 2 != len(self._bs.kLineSegs):
            self.print_warn("kpath string and extracted data are inconsistent. Skip")
            return
        locs = []
        labels = []
        # draw first and last symbol
        for i in [0, -1]:
            ksym = ksyms[i]
            s = KSYMBOL_LATEX.get(ksym, ksym)
            # coord = abs(i)
            coord = self._xs[i][i]
            # self._axes.annotate(s, xy=(coord, 0), xycoords="axes fraction", ha="center")
            locs.append(coord)
            labels.append(s)
        # draw intermediate points
        for i, x in enumerate(self._xs):
            if i == len(self._xs) - 1:
                break
            ksymLeft = ksyms[2*i+1]
            ksymRight = ksyms[2*i+2]
            if ksymLeft == ksymRight:
                s = KSYMBOL_LATEX.get(ksymLeft, ksymLeft)
            else:
                s = KSYMBOL_LATEX.get(ksymLeft, ksymLeft) + "|" + \
                            KSYMBOL_LATEX.get(ksymRight, ksymRight)
            # coord = xs[i][-1] / xs[-1][-1]
            coord = x[-1]
            # self._axes.annotate(s, xy=(coord, 0), xycoords="axes fraction", ha="center")
            locs.append(coord)
            labels.append(s)
        self._axes.set_xticks(locs)
        self._axes.set_xticklabels(labels)
        self._drawnKsym = True
    
    def draw_proj(self, atom, proj, *bands, **kwargs):
        '''
        '''
        bs = self._bs
        if not bs.hasProjection:
            self.print_warn("No wave projection is available. Skip.")
            return
        raise NotImplementedError
    
    def show(self):
        '''
        '''
        import matplotlib.pyplot as plt
        # hide the xticks if the kpoints symbols are not drawn
        if not self._drawnKsym:
            self._axes.get_xaxis().set_ticks([])
        plt.show()

    def export(self, *args, **kwargs):
        '''Wrapper to pyplot.savefig

        Args:
            args, kwargs: arguments parsed to savefig
        '''
        import matplotlib.pyplot as plt
        plt.savefig(*args, **kwargs)
