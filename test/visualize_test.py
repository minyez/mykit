#!/usr/bin/env python3
# coding=utf-8

import unittest as ut

from mykit.core.bandstructure import random_band_structure
from mykit.core.dos import Dos
from mykit.visualizer import BSVisualizer, init


class test_band_structure_visualizer(ut.TestCase):

    bs = random_band_structure(2, 10, 20, 4, 9, hasProjection=False)
    bsProj = random_band_structure(1, 10, 20, 4, 9, hasProjection=True)
    
    def test_draw(self):
        bsv = BSVisualizer(self.bs, ispin=0)
        bsv.draw('vbm')
        self.assertFalse(bsv.drawnSym)
        bsv.set_elim(-1,1)
        bsv.set_title("random band structure")
        # try a fake path
        bsv.mark_ksymbols("A-B")
        self.assertTrue(bsv.drawnSym)
        self.assertEqual(bsv.ispin, 0)
        bsv.ispin = 1
        self.assertEqual(bsv.ispin, 1)
        self.assertRaisesRegex(TypeError, r"ispin should be int", \
            bsv.__setattr__, "ispin", 2.0)
        self.assertRaisesRegex(ValueError, r"spin channel index overflow. nspins = 2", \
            bsv.__setattr__, "ispin", 3)
        self.assertRaisesRegex(TypeError, r"alignAtVbm should be bool", \
            bsv.__setattr__, "alignAtVbm", 1)
        
    def test_draw_proj(self):
        bsv = BSVisualizer(self.bs)
        bsvProj = BSVisualizer(self.bsProj)
        # raise if no projection is parsed
        self.assertRaisesRegex(AttributeError, \
            r"no projection data is available", \
            bsv.draw_proj, 'vbm', 'C', 's', 'vbm')
        bsvProj.draw_proj('vbm', (self.bsProj.atoms[0],), ('s',), 'vbm', label="test")

    def test_init(self):
        bsv = init(self.bs)
        self.assertIsInstance(bsv, BSVisualizer)
    

if __name__ == '__main__':
    ut.main()
