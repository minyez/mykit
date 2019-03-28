#!/usr/bin/env python3
# coding=utf-8

import unittest as ut

from mykit.core.bandstructure import _random_band_structure
from mykit.core.dos import Dos
from mykit.visualizer import BSVisualizer, init


class test_band_structure_visualizer(ut.TestCase):

    bs = _random_band_structure(1, 10, 20, 4, 9, hasProjection=False)
    bsProj = _random_band_structure(1, 10, 20, 4, 9, hasProjection=True)
    
    def test_draw(self):
        bsv = BSVisualizer(self.bs)
        bsv.draw('vbm')
        self.assertFalse(bsv.drawnSym)
        # try a fake path
        bsv.mark_ksymbols("A-B")
        # self.assertTrue(bsv.drawnSym)
    
    def test_draw_proj(self):
        bsv = BSVisualizer(self.bs)
        bsvProj = BSVisualizer(self.bsProj)
        # raise if no projection is parsed
        self.assertRaisesRegex(AttributeError, \
            r"no projection data is available", \
            bsv.draw_proj, 'vbm', 'C', 's', 'vbm')
        bsvProj.draw_proj('vbm', (self.bsProj.atoms[0],), ('s',), 'vbm')

    def test_init(self):
        bsv = init(self.bs)
        self.assertIsInstance(bsv, BSVisualizer)
    

if __name__ == '__main__':
    ut.main()
