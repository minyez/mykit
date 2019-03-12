# -*- coding: utf-8 -*-
'''this module defines some common used utilities'''


import os
import re


def get_dirpath(filePath):
    '''get the name of directory with filePath

    Args:
        filePath (str): the string of the path of file

    Returns:
        str: the path of parent directory, if filepath represents a file,
            otherwise the path of the directory
    '''

    _path = os.path.abspath(filePath)
    if os.path.isdir(_path):
        _path = _path + '/'
    return os.path.dirname(_path)


def trim_after(string, regex):
    '''Trim a string after the first match of regex.

    The matched pattern is trimed as well.
    
    Args:
        string (str): the string to trim
        regex (regex): the regex to match
    '''
    _search = re.search(regex, string)
    if _search != None:
        return string[:_search.start()]
    return string


def trim_before(string, regex):
    '''Trim a string from the beginning to the first match of regex.

    The matched pattern is not trimed.

    Args:
        string (str): the string to trim
        regex (regex): the regex to match
    '''
    _search = re.search(regex, string)
    if _search != None:
        return string[_search.start():]
    return string


def check_duplicates_in_tag_tuple(tagtuple):
    '''Check if there is duplicate in a tag tuple, case sensitive
    
    Args:
        tagTuple (tuple) : the tag tuple to check
    '''
    _dup = -1
    for _i, _k in enumerate(tagtuple):
        if _k in tagtuple[:_i]:
            _dup = _i
            break
    return _dup


def data_normalization(data, scale=1.0, normByPeak=True):
    '''Normalize the 1D data.

    Args:
        data (iterable): the container of 1D data
        normByPeak (bool) : when set True, the normalization factor will be
            the peak absolute value. Otherwise, the sum of absolute values 
            will be used as normalization factor.

    TODO:
        Generalize the normalization

    Returns:
        numpy array, the normalized data
    '''
    import numpy as np
    assert len(np.shape(data)) == 1
    assert isinstance(normByPeak, bool)
    
    _a = []
    try:
        _a = np.array(data, dtype="float64")
    except ValueError:
        raise ValueError("the data cannot be converted.")
    
    _sum = np.sum(np.abs(_a)) / scale
    _max = np.max(np.abs(_a)) / scale
    if normByPeak:
        return _a / _max
    return _a / _sum 


def find_data_extreme(data):
    '''Find the point at which the data reaches extrema

    TODO:
        Generalize to 2D and 3D coordinate

    Returns:
        dict, with two keys, "min" and "max".
        Either key has a 2-member tuple with its first the min/max value
        and second the coordinate where it reaches the extreme
    '''
    pass


def find_vol_dirs(path='.'):
    '''Find names of directories corresponding to calculation with lattice of different volumes
    
    The searching pattern is "V_x.xx" where x is 0-9
    '''
    pat = r'^V_\d.\d+'
    _dirs = []
    for _d in os.listdir(path):
        if re.match(pat, _d):
            _dirs.append(_d)
    def __sort_vol(dirstr):
        return float(dirstr.split('_')[1])
    _dirs = sorted(_dirs, key=__sort_vol)
    return _dirs
