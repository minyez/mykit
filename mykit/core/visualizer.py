# coding = utf-8

from mykit.core.bandstructure import BandStructure, init_band_visualizer
from mykit.core.dos import Dos


def init(*obj, plotter=None, **kwargs):
    '''Initialize a visualizer
    '''
    if len(obj) == 0:
        raise ValueError("should parse at least one object")
    if len(obj) == 1:
        o = obj[0]
        if isinstance(o, BandStructure):
            return init_band_visualizer(o, plotter=plotter, **kwargs)
        if isinstance(o, Dos):
            raise NotImplementedError
    if len(obj) == 2:
        pass
    raise NotImplementedError
