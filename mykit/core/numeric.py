# -*- coding: utf-8 -*-
'''Module that defines classes related to numerical issues
'''

from mykit.core.config import config

npPrec = config.get('numpyPrec')
del config

#pylint: disable=too-few-public-methods
class Prec:
    '''Metaclass to define the precision of numerical calculation and threshold
    '''
    _dtype = 'float'+str(npPrec).strip()
    _intType = 'int'+str(npPrec).strip()
    _symprec = 1.0E-5
    _thresEmp = 1.0E-3
    _thresOcc = 1.0 - _thresEmp
