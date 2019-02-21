# -*- coding: utf-8 -*-
'''
Module that defines class and utilities for band structure
'''

import warnings as wn
import numpy as np
from mykit.core.numeric import prec

class band_structure(prec):
    '''Band structure class.

    The energy bands should be parsed in a shape of (nspin,nkp,nbands).
    The dimension can be specified by corresponding arguments, which can also
    be set as 0 for automatic detection from the shape of ``bandE`` input.
    The projection of band wave function can also be parsed, with one extra
    
    Args:
      bandE (array-like) : the energy of all bands to be considered
      nspin (int) : 1 for non-sp calculation, 2 for sp calculation. 
      nband (int) : number of bands.
      nkp (int) : number of kpoints. 
      bandProj (array-like) : the projection of band wave function on atomic orbitals
      projName (list) : list of the atomic orbitals used
    '''
    def __init__(self, bandE, nspin=0, nbands=0, nkp=0, bandProj=None, projName=None):
        try:
            self.bandE = np.array(bandE, dtype=self._dtype)
            if not bandProj is None:
                self.bandProj = np.array(bandProj, dtype=self._dtype)
        except ValueError:
            raise ValueError("Non-array bandE or bandProj input")

        # check shape consistency
        wnStr = "Inconsistent {}: {} specified, {} in band shape"
        __eshape = np.shape(self.bandE)
        assert len(__eshape) == 3
        if nspin != 0:
            try:
                assert nspin == __eshape[0]
            except AssertionError:
                wn.warn(wnStr.format("nspin", nspin, __eshape[0]), wn.UserWarning)
        self.nspin = __eshape[0]

        if nbands != 0:
            try:
                assert nbands == __eshape[1]
            except AssertionError:
                wn.warn(wnStr.format("nbands", nbands, __eshape[1]), wn.UserWarning)
        self.nbands = __eshape[1]

        if nkp != 0:
            try:
                assert nkp == __eshape[2]
            except AssertionError:
                wn.warn(wnStr.format("nkp", nkp, __eshape[2]), wn.UserWarning)
        self.nkp = __eshape[2]

        if not bandProj is None:
            __pshape = np.shape(self.bandProj)
            assert len(__pshape) == 4
            assert __pshape[:3] == (self.nspin, self.nbands, self.nkp)
            assert __pshape[-1] == len(projName)





