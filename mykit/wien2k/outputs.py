# coding = utf-8

import os
from mykit.wien2k.utils import get_casename

class Vsp:
    """class for reading spherical part of potential with RMT in '.vsp' file

    Args:
        pathVsp (str): the file name of vsp file
    """

    def __init__(self, pathVsp=None):
        p = pathVsp
        if p is None:
            p = get_casename() + '.vsp'
        if not os.path.isfile(p):
            raise FileNotFoundError
        