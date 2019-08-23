# coding = utf-8
"""base class for unit conversion"""

from mykit.core.constants import ANG2AU, EV2HA, EV2RY, RY2HA


class UnitError(Exception):
    pass


class EnergyUnit:
    """Base class for controlling energy unit

    Args:
        eunit (str)
    """

    _defaultEU = 'ev'
    _validEU = ['ev', 'ry', 'au']
    _convEU = {
        ('ev', 'ry'): EV2RY,
        ('ev', 'au'): EV2HA,
        ('ry', 'au'): RY2HA,
    }

    def __init__(self, eunit=None):
        if eunit is None:
            self._eunit = self._defaultEU
        else:
            self._check_valid_eunit(eunit)
        self._eunit = eunit.lower()

    def _get_eunit_conversion(self, unitTo):
        self._check_valid_eunit(unitTo)
        tu = unitTo.lower()
        fu = self._eunit
        pair = (fu, tu)
        co = 1
        if pair in self._convEU:
            co = self._convEU[pair]
        elif pair[::-1] in self._convEU:
            co = 1.0 / self._convEU[pair[::-1]]
        return co
    
    def _check_valid_eunit(self, eunit):
        try:
            assert isinstance(eunit, str)
            u = eunit.lower()
            assert u in self._validEU
        except AssertionError:
            raise UnitError("allowed energy unit {}, {} parsed".format(
                self._validEU, eunit))


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

    def __init__(self, lunit=None):
        if lunit is None:
            self._lunit = self._defaultLU
        else:
            self._check_valid_lunit(lunit)
        self._lunit = lunit.lower()

    def _get_lunit_conversion(self, unitTo):
        self._check_valid_lunit(unitTo)
        tu = unitTo.lower()
        fu = self._lunit
        pair = (fu, tu)
        co = 1
        if pair in self._convLU:
            co = self._convLU[pair]
        elif pair[::-1] in self._convLU:
            co = 1.0 / self._convLU[pair[::-1]]
        return co

    def _check_valid_lunit(self, lunit):
        try:
            assert isinstance(lunit, str)
            u = lunit.lower()
            assert u in self._validLU
        except AssertionError:
            info = "allowed length unit {}, {} parsed".format(
                self._validLU, lunit)
            raise UnitError(info)
