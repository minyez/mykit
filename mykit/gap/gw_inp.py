# -*- coding: utf-8 -*
import re

from mykit.core.utils import trim_after


class GwInp():
    '''
    The changable parameters are saved in a dict, with the names of parameters as keys.
    The values are 6-member tuple, with members as

    1. Data type
    2. Pattern to match target block
    3. Default value
    4. Line in the block, 0 for one-line parameter
    5. Number of paramters in the block line, 0 for one-line parameter
    6. Index of parameter in the line, 0 for one-line parameter

    Args:
        path_gw_inp (str) : path to gw.inp
    '''
    available_params = {
        "pwm": (float, r"%BareCoul", 2.0, 1, 2, 0),
        "kmr": (float, r"%MixBasis", 0.75, 1, 1, 0),
        "barcevtol": (float, r"barcevtol", 0.1, 0, 0, 0),
        "barcevtol_soft": (float, r"barcevtol_soft", 0.1, 0, 0, 0),
        "MB_emax": (float, r"MB_emax", 20.0, 0, 0, 0),
        "lmbmax": (int, r"%MixBasis", 3, 2, 3, 0),
        "wftol": (float, r"%MixBasis", 1.0E-4, 2, 3, 1),
        "lblmax": (int, r"%MixBasis", 0, 2, 3, 2),
        "iop_core": (int, r"iop_core", 0, 0, 0, 0),
        "iop_fgrid": (int, r"%FreqGrid", 3, 1, 5, 0),
        "nomeg": (int, r"%FreqGrid", 16, 1, 5, 1),
        "omegmax": (float, r"%FreqGrid", 0.42, 1, 5, 2),
        "omegmin": (float, r"%FreqGrid", 0.00, 1, 5, 3),
        "emaxpol": (float, r"emaxpol", 1.0E10, 0, 0, 0),
        "emaxsc": (float, r"emaxsc", 1.0E10, 0, 0, 0),
        "n_ang_grid": (int, r"n_ang_grid", 6, 0, 0, 0),
        "lmax_q0": (int, r"lmax_q0", 6, 0, 0, 0),
        "iop_aniso": (int, r"iop_aniso", -1, 0, 0, 0),
        "iop_acfd": (int, r"%acfd", 1, 1, 1, 0),
        "exclude_cc_exx": (str, r"exclude_cc_exx", "F", 0, 0, 0), 
        }
    
    @classmethod
    def get_available_params(cls):
        """return the changable parameters"""
        return tuple(cls.available_params.keys())

    def __init__(self, path_gw_inp='gw.inp'):
        with open(path_gw_inp, 'r') as h:
            self._lines = [l for l in h.readlines()]
        self._params = {}
        self._locate_params()

    def _locate_params(self):
        '''locate line indices of parameters
        '''
        for i, line in enumerate(self._lines):
            l = trim_after(line, r'#').strip()
            if l == '':
                continue
            for k, v in self.available_params.items():
                if l.startswith(v[1]):
                    self._params[k] = i + v[3]
                    continue
    
    @property
    def params(self):
        '''dict, line index of parameter'''
        return self._params

    def __getitem__(self, key):
        return self.get_param(key)

    def get_param(self, key):
        '''Get the value of the parameter specified by key
        '''
        if key not in self.available_params:
            raise KeyError("%s is not available" % key)
        try:
            # extract from params
            l = self._params[key]
            t, _, _, _, n, i = self.available_params[key]
            if n == 0:
                # one-line parameter
                v = t(self._lines[l].split()[2])
            else:
                v = t(self._lines[l].split('l')[i])
        except KeyError:
            # return default value
            v = self.available_params[key][2]
        return v

    def modify_params(self, **kwargs):
        '''Change parameters

        Note:
            To modify parameters in block, the block should be present
            in the input file, otherwise it will be written as one-line
            parameter.
        '''
        gwlines = []
        linos = []
        params = []
        extras = []
        for k in kwargs:
            if k in self._params:
                params.append(k)
                linos.append(self._params[k])
            else:
                extras.append("%s = %s\n" % (k, kwargs[k]))
        for i, l in enumerate(self._lines):
            l = l.strip()
            if i in linos:
                k = params[linos.index(i)]
                pat = _get_pattern(self.available_params[k][-2])
                sub = _get_substr(*self.available_params[k][-2:], kwargs[k])
                l = re.sub(pat, sub, l.strip())
            gwlines.append(l+'\n')
        gwlines.extend(extras)
        return gwlines


_COMMENT_PAT = r"(#[\w \|\(\),\.-]*)?"


def _get_pattern(n):
    '''Return the pattern of n parameter block line'''
    if n > 0:
        s = [r"([\w \.-]+)", ] * n
        return r'^' + r'\|'.join(s) + _COMMENT_PAT + r'$'
    return r"^" + r"([\w \.-]+)=([\w \.-]+)" + _COMMENT_PAT + r"$"


def _get_substr(n, ind, value):
    '''Substitute parameter with value'''
    if n == 0:
        sub = "\\1 = " + str(value) + ' \\3'
    else:
        sublist = ["\\"+str(i+1) for i in range(n)]
        sublist[ind] = str(value)
        sub = ' | '.join(sublist) + ' \\' + str(n+1)
    return sub
