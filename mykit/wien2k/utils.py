# coding=utf-8

import os
from shutil import which
from fnmatch import fnmatch
from mykit.wien2k.constants import (DEFAULT_R0, DEFAULT_R0S, DEFAULT_RMT,
                                    DEFAULT_RMTS)
from mykit.core.utils import get_dirpath


def get_default_rmt_r0(elem):
    '''get the default RMT and R0 for element ``elem``
    '''
    rmt = DEFAULT_RMTS.get(elem, DEFAULT_RMT)
    r0 = DEFAULT_R0S.get(elem, DEFAULT_R0)
    return rmt, r0


def get_casename(w2kdir='.'):
    '''return the case name of a wien2k working directory
    
    It will first search for case.struct, if exists return the filename without extension 
    Otherwise the name of the directory will be returned
    '''
    abspath = os.path.abspath(w2kdir)
    if os.path.isdir(abspath):
        for filename in os.listdir(abspath):
            if fnmatch(filename, abspath + '/*.struct'):
                case = filename.split('/')[-1][:-8]
                return case

    return os.path.basename(get_dirpath(abspath))


def find_complex_file(casename, ext):
    '''Check if "{casename}.{ext}" or "{casename}.{ext}c" exists

    Args:
        casename (str) 
        ext (str) : the file extension without "c"
    '''
    _path = casename + '.' + ext
    if not os.path.isfile(_path):
        _path = _path + 'c'
    if not os.path.isfile(_path):
        raise FileNotFoundError("Neither {casename}.in1 nor {casename}.in1c is found.")
    return _path


def get_run_lapw():
    '''Return the absolute path of run_lapw executable

    Returns:
        str
    '''
    return which('run_lapw')
