# -*- coding: utf-8 -*-
'''Module that defines classes related to numerical issues
'''

from mykit.core.config import global_config

#pylint: disable=too-few-public-methods
class prec:
    '''Metaclass to define the precision of numerical calculation and threshold
    '''
    _dtype = global_config.get('numpyPrec')
    # _order = 'C'
    _symprec = 1.0E-5

    # def __init__(self, dtype=None, order=None, symprec=None):
    #     if dtype:
    #         self._dtype = dtype
    #     if order:
    #         self._order = order
    #     if symprec:
    #         self._symprec = symprec