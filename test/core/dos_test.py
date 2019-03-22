#!/usr/bin/env python3
# coding = utf-8

import unittest as ut

import numpy as np

from mykit.core.dos import Dos, DosError

goodEgrid = [-4,-3,-2,-1,0,1,2,3,4]
badEgrid = [-4,-3,-2,-1,0,1]
goodDos = [
    [1.0, 1.0],
    [2.0, 2.0],
    [3.0, 4.0],
    [1.0, 2.0],
    [0.0, 0.0],
    [0.0, 0.0],
    [0.0, 0.0],
    [1.0, 2.0],
    [2.0, 1.0],
]
badDos1 = [1.0,2.0,3.0,4.0],
badDos2 = [
    [3.0, 4.0],
    [1.0, 2.0],
    [0.0, 0.0],
    [0.0,]
]
badDos = (badDos1, badDos2)
nedos, nspins = np.shape(goodDos)


class test_dos_initialize(ut.TestCase):

    def test_raise_for_inconsistent_egrid_dos(self):
        for b in badDos:
            self.assertRaisesRegex(DosError, r"Inconsistent shape: *", \
                Dos, goodEgrid, b, 0.0)
        self.assertRaisesRegex(DosError, r"Inconsistent shape: *", \
            Dos, badEgrid, goodDos, 0.0)
    
    def test_properties(self):
        dos = Dos(goodEgrid, goodDos, efermi=0.0)
        self.assertEqual(dos.nedos, nedos)
        self.assertEqual(dos.nspins, nspins)
        



if __name__ == '__main__':
    ut.main()
