# coding = utf-8
"""base class for unit conversion"""

from mykit.core.constants import ANG2AU, EV2HA, EV2RY, RY2HA


class UnitError(Exception):
    pass

class _unit:
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
    # add reverse conversion
    pairs = list(_convEU.keys())
    for p in pairs:
        _convEU[p[::-1]] = 1.0 / _convEU[p]

    def __init__(self, eunit=None):
        if eunit is None:
            self._eunit = self._defaultEU
        else:
            self._check_valid_eunit(eunit)
        self._eunit = eunit.lower()

    def _get_eunit_conversion(self, eunit):
        self._check_valid_eunit(eunit)
        tu = eunit.lower()
        fu = self._eunit
        pair = (fu, tu)
        co = 1
        if pair in self._convEU:
            co = self._convEU[pair]
        return co
    
    def _check_valid_eunit(self, eunit):
        try:
            assert isinstance(eunit, str)
            u = eunit.lower()
            assert u in self._validEU
        except AssertionError:
            raise UnitError("allowed energy unit {}, {} parsed".format(
                self._validEU, eunit))

    @classmethod
    def convert_eunit(cls, value, unit_from, unit_to):
        """convert `value` in unit `u_from` to unit `u_to`

        Args:
            value (number or ndarray)
            unit_from (str)
            unit_to (str)
        """
        try:
            assert unit_from in cls._validEU
            assert unit_to in cls._validEU
        except AssertionError:
            raise UnitError("unit not supported")
    
        t = (unit_from ,unit_to)
        if t in cls._convEU:
            factor = cls._convEU[t]
        else:
            raise UnitError("Conversion relation not set for %s and %s" % (unit_from, unit_to))
        return value * factor

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
    # add reverse conversion
    pairs = list(_convLU.keys())
    for p in pairs:
        _convLU[p[::-1]] = 1.0 / _convLU[p]

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

    @classmethod
    def convert_lunit(cls, value, unit_from, unit_to):
        """convert `value` in unit `u_from` to unit `u_to`

        Args:
            value (float)
            unit_from (str)
            unit_to (str)
        """
        try:
            assert unit_from in cls._validLU
            assert unit_to in cls._validLU
        except AssertionError:
            raise UnitError("unit not supported")
    
        t = (unit_from ,unit_to)
        if t in cls._convLU:
            factor = cls._convLU[t]
        else:
            raise UnitError("Conversion relation not set for %s and %s" % (unit_from, unit_to))
        return value * factor
