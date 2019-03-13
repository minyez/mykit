# coding = utf-8

import os
import shutil
import subprocess as sp
from copy import deepcopy
from fnmatch import fnmatch

from mykit.core.config import global_config
from mykit.core.log import verbose


class PotcarError(Exception):
    pass


class potcar_search(verbose):
    '''Class for locating POTCAR.

    Args:
        names (str): the names of elements to search
        usegw (bool): if set True, only POTCARs with GW suffix will be searched or exported
    '''

    def __init__(self, *names, usegw=False):
        assert isinstance(usegw, bool)
        config = global_config()
        _homePawPbe, _homePawLda = config.get("vaspPawPbe", "vaspPawLda")
        del config
        self._homePaw = {"PBE": _homePawPbe, "LDA": _homePawLda}
        try:
            assert len(names) > 0
        except AssertionError:
            raise PotcarError("should have at least one element")
        
        self._usegw = usegw
        self._names = names

    @property
    def names(self):
        _names = list(deepcopy(self._names))
        for i, name in enumerate(_names):
            if self._usegw:
                if not name.endswith("_GW"):
                    _names[i] = name + "_GW"
        return tuple(_names)
    @names.setter
    def names(self, value):
        try:
            assert isinstance(value, (str, list, tuple))
            if isinstance(value, (list, tuple)):
                for _i, v in enumerate(value):
                    assert isinstance(v, str)
        except AssertionError:
            raise PotcarError("invalid name, should be (list/tuple of) string : {}".format(value))
        if isinstance(value, str):
            self._names = (value,)
        else:
            self._names = (*value,)

    def _get_xc_home(self, xc):
        '''Get the home directory of xc type PAW

        Args:
            xc (str): name of xc, case insensitive.
        '''
        if xc is None:
            _xc = "PBE"
        else:
            try:
                assert isinstance(xc, str)
            except AssertionError:
                raise PotcarError("XC should be string type")
            _xc = xc.upper()
        try:
            assert _xc in self._homePaw
        except AssertionError:
            raise PotcarError("XC type not supported: {}".format(_xc))
        _home = self._homePaw[_xc]
        if _home is None or not os.path.isdir(_home):
            raise PotcarError("{} PAW home directory does not exist: {}".format(_xc, _home))
        return _home

    def search(self, xc="PBE"):
        '''Search the PAW xc directory for available POTCARs for specified element names

        Args:
            xc (str): name of xc, case insensitive.
        
        Returns:
            str, the path of the searching directory
            dict, items are name-dict pair, with each item in dict is POTCAR:(ENMIN, ENMAX)
        '''
        _home = self._get_xc_home(xc)
        _dict = {}
        for name in self._names:
            _d = {}
            for i in os.listdir(_home):
                if fnmatch(i, name) or fnmatch(i, name+"_*"):
                    if not i.endswith("_GW") and self._usegw:
                        continue
                    _path = os.path.join(_home, i, 'POTCAR')
                    enmax = get_enmax(_path)
                    enmin = get_enmin(_path)
                    _d.update({i: (enmin, enmax)})
            _dict.update({name: _d})
        return _home, _dict
        
    def export(self, xc="PBE", pathPotcar='POTCAR'):
        '''Concentate POTCARs of element names to ``pathPotcar``

        Args:
            xc (str): name of xc, case insensitive.
            pathPotcar (str): the path to output the concentated POTCAR.

        Returns:
            str, the searching directory
        '''
        _home = self._get_xc_home(xc)
        # ! The following exception is not unittest documented
        _files = []
        for _i, name in enumerate(self.names):
            path = os.path.join(_home, name, 'POTCAR')
            if not os.path.isfile(path):
                raise PotcarError("POTCAR not found: {}".format(path))
            _files.append(path)
        with open(pathPotcar, 'w') as h_o:
            for f in _files:
                with open(f, 'r') as h_i:
                    for line in h_i:
                        h_o.write(line)
        return _home


def get_enmax(pathPotcar):
    '''Get ENMAX of a POTCAR
    '''
    return float(sp.check_output(['grep', 'ENMAX', pathPotcar]).split()[2][:-1])


def get_enmin(pathPotcar):
    '''Get ENMIN of a POTCAR
    '''
    return float(sp.check_output(['grep', 'ENMIN', pathPotcar]).split()[5])
