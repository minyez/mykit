# -*- coding: utf-8 -*-
'''Module that defines classes related to numerical issues
'''

from mykit.core.config import global_config

config = global_config()
npPrec = config.get('numpyPrec')
del config, global_config

#pylint: disable=too-few-public-methods
class Prec:
    '''Metaclass to define the precision of numerical calculation and threshold
    '''
    _dtype = 'float'+str(npPrec).strip()
    _symprec = 1.0E-5
