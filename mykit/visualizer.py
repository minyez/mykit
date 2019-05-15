# coding = utf-8

from copy import deepcopy

import numpy as np

from mykit.core.bandstructure import BandStructure
from mykit.core.dos import Dos
from mykit.core.kmesh import KSYMBOL_LATEX, KmeshError, kpath_decoder
from mykit.core.log import Verbose


def init(*objs, **kwargs):
    '''Convenient function to initialize a visualizer

    Args:
        objs: the objects to parse to the visualizer.
        kwargs: keyword arguments to parse to the visualizer.
    '''
    if len(objs) == 0:
        raise ValueError("should parse at least one object")
    if len(objs) == 1:
        o = objs[0]
        if isinstance(o, BandStructure):
            return BSVisualizer(o, **kwargs)
        if isinstance(o, Dos):
            raise NotImplementedError
    if len(objs) == 2:
        pass
    raise NotImplementedError


class BSVisualizer(Verbose):
    '''The class for drawing band structure with the help of matplotlib.pyplot.

    Currently only support drawing one bandstructure, i.e. 1 Axes in Fig.
    Can only plot one spin channel.

    Args:
        bs (``BandStructure``): the band structure to plot
        ispin (int): the spin channel to show.
        dos (`Dos` object or dict): density of states to draw with the band structure
            if a dict is parsed, it will be used as the keyword arguments for ``get_dos``
            method of ``BandStructure``
        projStyle ('dot', 'stripe'): the style to draw the wave projections
        kwargs: keyword arguments to parse to initialize with `pyplot.subplots`

    TODO:
        simultaneous drawing with xmgrace object
    '''
    _SUPPORT_PROJ_STYLE = ['dot', 'stripe']

    def __init__(self, bs, align_vbm=False, dos=None, proj_style='dot', ispin=0,
                 **kwargs):
        assert isinstance(bs, BandStructure)
        # argument type check
        if not bs.isKpath:
            self.print_warn(
                "k-vectors of Band structure do not form a kpath. Plotting anyway...")
        try:
            import matplotlib.pyplot as plt
        except ModuleNotFoundError:
            raise ValueError("Matplotlib is required to use pyplot plotter")
        self._drawDos = False
        if dos is not None:
            self._drawDos = True
            if isinstance(dos, dict):
                self._dos = bs.get_dos(**dos)
            elif isinstance(dos, Dos):
                _check_bs_dos_consistency(bs, dos)
                self._dos = deepcopy(dos)
            else:
                raise TypeError(
                    "dos should be a number or Dos objest")

        if proj_style not in self._SUPPORT_PROJ_STYLE:
            raise ValueError("projStyle {} is not supported. {}".format(
                proj_style, self._SUPPORT_PROJ_STYLE))

        self._bs = deepcopy(bs)
        self._xs = bs._generate_kpath_x()
        self._efermi = bs.efermi
        self._nspins = bs.nspins
        self.ispin = ispin
        self.alignAtVbm = align_vbm
        self._drawnKsym = False
        self._projStyle = proj_style
        self._useLabel = False
        # initialize the Figure and Axes
        if self._drawDos:
            kwargs.pop("nrows", None)
            kwargs.pop("ncols", None)
            kwargs.pop("sharey", None)
            self._fig, (self._axbs, self._axdos) = plt.subplots(
                1, 2, gridspec_kw={"width_ratios": [2, 1]}, sharey=True, **kwargs)
            self._fig.subplots_adjust(wspace=0)
        else:
            self._fig, self._axbs = plt.subplots(**kwargs)
        # set band structre axis
        self._axbs.set_xlim([self._xs[0][0], self._xs[-1][-1]])
        # set x tick length to zero to make them invisible
        self._axbs.tick_params(axis=u'x', which=u'both', length=0)
        self._axbs.set_ylabel("Energy ({})".format(bs.unit))
        # set dos axis
        if self._drawDos:
            # hide marks and tick of x axis
            self._axdos.get_xaxis().set_ticks([])
            self._axdos.set_xlim([0.0, np.max(self._dos.dos)*1.1])
            self._draw_total_dos()

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
            raise ValueError(
                "spin channel index overflow. nspins = %d" % self._nspins)
        self._ispin = newValue

    def set_title(self, title, **kwargs):
        '''Set the title of figure

        Wrapper of pyplot.Axes.set_title
        '''
        self._axbs.set_title(title, **kwargs)

    def set_elim(self, bottom, top, **kwargs):
        '''Set the energy limit to show

        Wrapper of Axes.set_ylim

        Args:
            bottom (float): the lower bound of the energy range
            top (float): the upper bound of the energy range
            kwargs: the keyword arguments to parse to Axes.set_ylim
        '''
        self._axbs.set_ylim(bottom=bottom, top=top, **kwargs)

    def _draw_total_dos(self):
        '''draw the total dos
        '''
        dos = self._dos
        edos = dos.edos - dos.efermi * int(self.alignAtVbm)
                    # self._bs.vbmPerSpin[self.ispin] * int(self.alignAtVbm)
        self._axdos.plot(dos.dos[:, self.ispin], edos, color="k")
        # draw Fermi level
        if self.alignAtVbm:
            self._axdos.axhline(0.0, color="k", linestyle='dashed')
        else:
            self._axdos.axhline(dos.efermi, color="k", linestyle='dashed')


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
        b = bs.get_band_indices(*bands)
        for i, (stk, edk) in enumerate(bs.kLineSegs):
            for ib in b:
                eigen = bs.eigen[self.ispin, stk:edk+1, ib] - \
                    bs.vbmPerSpin[self.ispin] * int(self.alignAtVbm)
                    # bs.vbmPerSpin[self.ispin] * int(self.alignAtVbm)
                self._axbs.plot(xs[i], eigen, **kwargs)
            # draw vertical line to separate the different line segments
            if i != len(bs.kLineSegs) - 1:
                self._axbs.axvline(xs[i][-1], color="k")
        # draw Fermi level
        if self.alignAtVbm:
            self._axbs.axhline(0.0, color="k", linestyle='dashed')
        else:
            self._axbs.axhline(bs.efermi, color="k", linestyle='dashed')

    def mark_ksymbols(self, kpath):
        '''Mark the kpoints symbols on the plot.

        Args
            kpath (str): the kpath string.
                The string will be decoded by `kpath_decoder` to a list of kpoints symbols.
                followed by a consistency check.
                If the check fail, a warning will be prompted and no symbols plotted
        '''
        try:
            ksyms = kpath_decoder(kpath)
        except KmeshError:
            return
        if len(ksyms) / 2 != len(self._bs.kLineSegs):
            self.print_warn(
                "kpath string and extracted data are inconsistent. Skip")
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
        self._axbs.set_xticks(locs)
        self._axbs.set_xticklabels(labels)
        self._drawnKsym = True

    def draw_proj(self, atom, proj, *bands, **kwargs):
        '''draw the wave projection on the band structure diagram

        Args:
            atom (int, str or Iterable):
            proj (int, str or Iterable):
            bands (int): the indices of bands to draw the projection
            kwargs: keyword argument to parse to pyplot.Axes.scatter or pyplot.Axes.fill_between
                depending on the projStyle set at initialization.
        '''
        bs = self._bs
        xs = self._xs
        amplifier_dot = 120.0
        # use triple band gap as multiplier for stripe mode
        amplifier_stripe = bs.fundGap[self.ispin] * 3
        if not bs.hasProjection:
            raise AttributeError("no projection data is available")
        # get projection data
        pWave = bs.sum_atom_proj_comp(atom, proj, fail_one=False)
        binds = bs.get_band_indices(*bands)
        if 'color' not in kwargs:
            kwargs['color'] = 'k'
        if 's' in kwargs:
            kwargs.pop('s')
        if 'label' in kwargs:
            self._useLabel = True
        # draw DOS
        if self._drawDos:
            dos = self._dos
            edos = dos.edos - dos.efermi * int(self.alignAtVbm)
            pDos = dos.sum_atom_proj_comp(atom, proj, fail_one=False)
            self._axdos.plot(pDos[:, self.ispin], edos, **kwargs)
        # draw band
        for i, (stk, edk) in enumerate(bs.kLineSegs):
            for _j, bi in enumerate(binds):
                eigen = bs.eigen[self.ispin, stk:edk+1, bi] - \
                    bs.vbmPerSpin[self.ispin] * int(self.alignAtVbm)
                if self._projStyle == 'dot':
                    s = pWave[self.ispin, stk:edk+1, bi] * amplifier_dot
                    self._axbs.scatter(xs[i], eigen, s=s, **kwargs)
                if self._projStyle == 'stripe':
                    self._axbs.fill_between(xs[i], eigen,
                                            eigen - pWave[self.ispin, stk:edk+1, bi] * amplifier_stripe, **kwargs)
                # pop the label keyword such that label is only added for once
                if 'label' in kwargs:
                    kwargs.pop('label')

    def show(self, legend=True, tight_layout=False):
        '''Preview the band structure diagram. 

        A wrapper for pyplot.legend and pyplot.show
        '''
        import matplotlib.pyplot as plt
        # # hide the xticks if the kpoints symbols are not drawn
        # if not self._drawnKsym:
        #     self._axes.get_xaxis().set_ticks([])
        if self._useLabel and legend:
            plt.legend()
        if tight_layout:
            plt.tight_layout()
        plt.show()

    def export(self, *args, **kwargs):
        '''Wrapper to pyplot.savefig

        TODO:
            export agr file

        Args:
            args, kwargs: arguments parsed to savefig
        '''
        import matplotlib.pyplot as plt
        plt.savefig(*args, **kwargs)


def _check_bs_dos_consistency(bs, dos):
    '''Check the consistency of BandStructure and Dos objects

    Spin, unit, projections (projectors and atoms)
    '''

    assert isinstance(bs, BandStructure)
    assert isinstance(dos, Dos)
    assert bs.nspins == dos.nspins
    assert bs.unit == dos.unit
    assert bs.hasProjection == dos.hasProjection

    if bs.hasProjection:
        assert bs.natoms == dos.natoms
        assert bs.nprojs == dos.nprojs
        for i, a in enumerate(bs.atoms):
            assert a == dos.atoms[i]
        for i, p in enumerate(bs.projs):
            assert p == dos.projs[i]
