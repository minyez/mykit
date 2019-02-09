# -*- coding: utf-8 -*-
'''this module defines some common used utilities'''


import os

def get_dirpath(filepath):
    '''get the name of directory with filepath

    Args:
        filepath (str): the string of the path of file

    Returns:
        str: the path of parent directory if filepath represents a file,
            otherwise the path of the directory
    '''

    __path = os.path.abspath(filepath)
    if os.path.isdir(__path):
        __path = __path + '/'
    return os.path.dirname(__path)
