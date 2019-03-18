# -*- coding: utf-8 -*-
'''include class for KPOINTS
'''

import os
from collections.abc import Iterable

from mykit.core.kmesh import _check_valid_kpath_dict, kmesh_control
from mykit.core.log import verbose


class KpointsError(Exception):
    pass

# ! kmesh_control instance is used as an attribute,
# ! since vasp does not have any tag for it
class kpoints(verbose):
    '''class for manipulation KPOINTS, IBZKPT files

    At least one of kdiv, kdense, kpath, kpoints
    should be specified.

    Note:
        When auto mode check is triggered (kmode=None),
        the priority to determine the mode is:
            kpoints > kpath > kdiv > kdense

    Args:
        comment (str): the comment for KPOINTS
        kdiv (3-member list): the number of grid on each axis
        kdense (integer): the density of k grid point
            If set to non-zero value and kmode is "G", "M" or "A", the kmode will automatically be set as "A".
            If kmode is set as "L" and kpath specified, it is the number of kpoints along two special points.
        kshift (3-member list): shift off center
        kmode ("G", "M", "A", "L"): the mode of KPOINTS
            None: use explicit kpoints format
            "G": Gamma-centered
            "M": Monkhorse-Pack
            "A": fully automatic, i.e. the grid will be decided by kdense
            "L": line mode for bandstructure. kpath should be parse as well. Reciprocal vector unit is always used.
        kpath (dict): the ends of each line segment of a kpath. "symbols" mark the symbols of ends and "coordinates"
            give their corresponding coordinates. For example:
            ```
                kpath = {
                    "symbols": ["GM","L","L","GM],
                    "coordinates": [[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0]]
                }
            ```
            Both keys words should have a value of list with even length.
            If specified, the mode will be changed to "L".
        kpoints (nx4 array): the explicit kpoints to parse. the fourth column is the kpoint weight.
            If specified, the explict mode will be triggered.
    '''

    def __init__(self, comment=None, kmode=None, kdiv=None, kdense=0, kshift=None, kpath=None, kpoints=None):
        try:
            assert kdense >= 0
        except AssertionError:
            raise KpointsError("kdense must not be negative")

        self.comment = comment
        # self._relatePoscar = None
        _kmode = kmode
        if isinstance(_kmode, str):
            _kmode = _kmode.upper()
        try:
            assert _kmode in [None, "G", "M", "L", "A"]
        except AssertionError:
            raise KpointsError("Unknown KPOINTS mode: {}".format(_kmode))

        # Automatic mode check when kmode is not specified
        if _kmode == None:
            if kpoints != None:
                pass
            elif kpath != None:
                _check_valid_kpath_dict(kpath)
                _kmode = "L"
            elif kdiv != None:
                _kmode = "G"
            elif kdense > 0:
                # self.print_warn("Length(kdense) set for non-line mode. Switch to \"A\" mode.", level=1)
                _kmode = "A"
            else:
                raise KpointsError("Enter at least one of kdiv, kdense, kpath and kpoints")
        # Consistency check for each mode
        if _kmode == "L":
            if kdense == 0:
                raise KpointsError("kdense > 0 for line mode")
            if kpath == None:
                raise KpointsError("kpath should be specified for line mode")
        if _kmode == "A" and kdense == 0:
            raise KpointsError("Fully automatic mode needs positive kdense (e.g. 30)")
        if _kmode in ["G", "M"] and kdiv == None:
            raise KpointsError("kdiv should be specified for Gamma or MP mode")

        self._control = kmesh_control("vasp", kdiv=kdiv, kdense=kdense, \
            kshift=kshift, kpath=kpath, kpoints=kpoints, kmode=_kmode)

    @property
    def mode(self):
        return self._control.kmode

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        _strs = []
        if not self.comment:
            # TODO: automatic add comment about kpath in line mode. Or determine it in script
            _strs += ["KPOINTS generated by mykit",]
        else:
            _strs += [self.comment,]
        _explictKs = self._control.tag_vals("vasp", "kpoints")[0]
        if _explictKs != None:
            _strs += [len(_explictKs), "Reciprocal",]
            for _i, kp in enumerate(_explictKs):
                try:
                    _strs.append("%20.14f %20.14f %20.14f %14d" % kp)
                except TypeError:
                    raise KpointsError("Bad explicit kpoints format: {}, index {}".format(kp, _i))
        else:
            _modeDict = {"G": "Gamma", "M": "Monkhorse-Pack", "L": "Line", "A": "Auto"}
            _kdense = self._control.tag_vals("vasp", "kdense")[0]
            _mode = self._control.tag_vals("vasp", "kmode")[0]
            if _mode == "L":
                kd = int(_kdense)
            else:
                kd = 0
            _strs += [str(kd), _modeDict[_mode],]
            if _mode in ["G", "M"]:
                _kg, _ks = self._control.tag_vals("vasp", "kdiv", "kshift")
                try:
                    _strs.append("{:2d} {:2d} {:2d}".format(*_kg))
                except (TypeError, IndexError):
                    raise KpointsError("Bad kdiv format for G/M mode: {}".format(_kg))
                if _ks == None:
                    _ks = [0, 0, 0]
                try:
                    _strs.append("{:2d} {:2d} {:2d}".format(*_ks))
                except (TypeError, IndexError):
                    raise KpointsError("Bad kshift format for G/M mode: {}".format(_ks))
            elif _mode == "A":
                _strs.append(str(int(_kdense)))
            elif _mode == "L":
                _strs.append("Reciprocal")
                _kpath = self._control.tag_vals("vasp", "kpath")[0]
                try:
                    assert isinstance(_kpath, dict) and "symbols" in _kpath and "coordinates" in _kpath
                except AssertionError:
                    raise KpointsError("Invalid kpath as a dict")
                symbols = _kpath["symbols"]
                coords = tuple(_kpath["coordinates"])
                if len(symbols) != len(coords):
                    raise KpointsError("Inconsistent length of symbols and coordiantes")
                if len(symbols)%2 != 0:
                    raise KpointsError("Odd length found for symbols/coordiantes, require even.")
                nSeg = int(len(symbols)/2)
                for i in range(nSeg):
                    st = 2*i
                    ed = 2*i + 1
                    _strs.extend(["{:9.6f} {:9.6f} {:9.6f} {} #{}".format(*coords[st], 1, symbols[st]), \
                        "{:9.6f} {:9.6f} {:9.6f} {} #{}".format(*coords[ed], 1, symbols[ed]), \
                        ''])
            else:
                raise KpointsError("Unknown KPOINTS mode: {}".format(_mode))
        return '\n'.join(_strs)
        
    def write(self, pathKpoints="KPOINTS", backup=False, suffix="_bak"):
        '''Write KPOINTS to path
        '''
        _name = pathKpoints
        try: 
            assert not os.path.isdir(_name)
        except AssertionError:
            raise KpointsError("The path to write KPOINTS is a directory.")
        if os.path.isfile(_name) and backup:
            _bakname = _name + suffix.strip()
            os.rename(_name, _bakname)
        with open(_name, 'w') as f:
            print(self.__str__(), file=f)

    @classmethod
    def read_from_kpoints(cls, pathKpoints="KPOINTS"):
        pass

    @classmethod
    def read_from_ibzkpt(cls, pathKpoints="IBZKPT"):
        pass

    # * Factory methods to write kpath of particular crystal
    # * Name with kpath_(bravis)
