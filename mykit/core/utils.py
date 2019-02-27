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

def trim_comment(string, commMark):
    '''Trim a comment string from the end
    
    Args:
        string (str): the string to trim
        commMark (regex): the comment mark to trim from
    '''
    _search = re.search(commMark, string)
    if _search != None:
        return string[:_search.start()]
    else:
        return string


def check_duplicates_in_tag_tuple(tuple):
    '''Check if there is duplicate in a tag tuple, case sensitive
    
    Args:
        tagTuple (tuple) : the tag tuple to check
    '''
    _dup = -1
    for _i, _k in enumerate(tuple):
        if _k in tuple[:_i]:
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