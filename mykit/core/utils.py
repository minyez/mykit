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
