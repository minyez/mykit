# coding = utf-8

import numpy as np

from mykit.core.constants import EV2HA, EV2RY, RY2HA


class EnergyUnit:
    '''Base class for controlling energy unit

    Args:
        eunit (str)
    '''

    _defaultEU = 'ev'
    _validEU = ['ev', 'ry', 'au']
    _convEU = {
        ('ev', 'ry'): EV2RY,
        ('ev', 'au'): EV2HA,
        ('ry', 'au'): RY2HA,
    }

    def __init__(self, eunit):
        try:
            u = eunit.lower()
            assert u in self._validEU
        except (AttributeError, AssertionError):
            u = self._defaultEU
        self._eunit = u

    def _get_eunit_conversion(self, unitTo):    
        tu = unitTo.lower()
        fu = self._eunit
        pair = (fu, tu)
        co = 1
        if pair in self._convEU:
            co = self._convEU[pair]
        elif pair[::-1] in self._convEU:
            co = 1.0 / self._convEU[pair[::-1]]
        return co
