# coding = utf-8

import numpy as np

from mykit.core.constants import ANG2AU, EV2HA, EV2RY, RY2HA


class UnitError(Exception):
    pass


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
        try:
            assert isinstance(unitTo, str)
            assert unitTo.lower() in self._validEU
        except AssertionError:
            raise UnitError("allowed energy unit {}, {} parsed".format(self._validEU, unitTo))
        tu = unitTo.lower()
        fu = self._eunit
        pair = (fu, tu)
        co = 1
        if pair in self._convEU:
            co = self._convEU[pair]
        elif pair[::-1] in self._convEU:
            co = 1.0 / self._convEU[pair[::-1]]
        return co


class LengthUnit:
    '''Base class for controlling length unit

    Args:
        lunit (str)
    '''

    _defaultLU = 'ang'
    _validLU = ['ang', 'au']
    _convLU = {
        ('ang', 'au'): ANG2AU,
    }

    def __init__(self, lunit):
        try:
            u = lunit.lower()
            assert u in self._validLU
        except (AttributeError, AssertionError):
            u = self._defaultLU
        self._lunit = u

    def _get_lunit_conversion(self, unitTo):
        try:
            assert isinstance(unitTo, str)
            assert unitTo.lower() in self._validLU
        except AssertionError:
            raise UnitError("allowed length unit {}, {} parsed".format(self._validLU, unitTo))
        tu = unitTo.lower()
        fu = self._lunit
        pair = (fu, tu)
        co = 1
        if pair in self._convLU:
            co = self._convLU[pair]
        elif pair[::-1] in self._convLU:
            co = 1.0 / self._convLU[pair[::-1]]
        return co
