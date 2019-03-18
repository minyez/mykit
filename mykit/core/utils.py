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


def common_ss_conv(string, i, conv2, sep=None):
    '''
    Split the string and convert a single substring to a specified type.

    Args:
        string (str): the string from which to convert value
        i (int): the substring index in the list to be converted after splitting by sep
        conv2: the type to which the substring will be converted
        sep (regex): the separators used to split the string.
    '''

    str_tmp = string.strip()
    if sep is not None:
        str_list = re.split(r'[%s]'%sep, str_tmp)
    #    print(str_list)
    else:
        str_list = str_tmp.split()

    # try:
        # return conv2(str_list[i])
    # except ValueError:
    #     return conv2(float(str_list[i]))
    return conv2(str_list[i])


def get_first_last_line(filePath):
    '''Return the first and the last non-empty line of file

    The existence of filePath should be check first.

    Args:
        filePath (str): the path of the file
    '''
    from sys import stdin
    with open(filePath, "rb") as f:
        first = f.readline()        # Read the first line.
        f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
        while f.read(1) != b"\n":   # Until EOL is found...
            f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
        last = f.readline()         # Read last line.
    # convert to string by stdin encoding
    return str(first, stdin.encoding).strip(), str(last, stdin.encoding).strip()
