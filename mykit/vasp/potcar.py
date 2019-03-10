# coding = utf-8

import os
import shutil

from mykit.core.config import global_config
from mykit.core.log import verbose


class PotcarError(Exception):
    pass


class potcar_search(verbose):
    '''Class for locating POTCAR.

    Args:
        names (str) : the name of POTCARs to use. GW suffix can be added
        xc ("PBE", "LDA") : the xc type of POTCAR, , case-insensitive
        usegw (bool): if set True, the GW suffix is automatically added for all elements
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
        
        self._names = list(names)
        for i, name in enumerate(self._names):
            if usegw:
                if not name.endswith("_GW"):
                    self._names[i] = name + "_GW"

    @property
    def names(self):
        return self._names
    
    def export(self, xc="PBE", pathPotcar='POTCAR'):
        '''Concentate POTCAR of element names to ``pathPotcar``

        Args:
            pathPotcar (str): the path to output the concentated POTCAR.
        '''
        _xc = xc.upper()
        try:
            assert _xc in ["PBE", "LDA"]
        except AssertionError:
            raise PotcarError("Support PBE(PE) and LDA(CA) only")

        _home = self._homePaw[_xc]
        if _home is None or not os.path.isdir(_home):
            raise PotcarError("{} PAW home directory does not exist: {}".format(_xc, _home))

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
