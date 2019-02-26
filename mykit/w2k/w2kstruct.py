# coding = utf-8
'''module that defines class and functions that manipulates WIEN2k main input file struct
'''
from mykit.core.lattice import lattice


class structError(Exception):
    pass


class w2kstruct(lattice):
    '''Class for WIEN2k main input file
    '''
    casename = ''

    def __init__(self, cell, atoms, pos, **kwargs):
        super(w2kstruct, self).__init__(cell, atoms, pos, **kwargs)
