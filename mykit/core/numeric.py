# -*- coding: utf-8 -*-
'''Module that defines classes related to numerical issues
'''

# TODO: import settings from a metadata file
#pylint: disable=too-few-public-methods
class prec:
    '''Metaclass to define the precision of numerical calculation and threshold
    '''
    _dtype = 'float32'
    _order = 'C'
    _symprec = 1.0E-5

