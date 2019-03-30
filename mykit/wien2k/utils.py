# coding=utf-8

from mykit.wien2k.constants import (DEFAULT_R0, DEFAULT_R0S, DEFAULT_RMT,
                                    DEFAULT_RMTS)


def get_default_rmt_r0(elem):
    '''get the default RMT and R0 for element ``elem``
    '''
    rmt = DEFAULT_RMTS.get(elem, DEFAULT_RMT)
    r0 = DEFAULT_R0S.get(elem, DEFAULT_R0)
    return rmt, r0
