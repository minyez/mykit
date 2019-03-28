#!/usr/bin/env python3
# coding = utf-8

import unittest as ut

import numpy as np

from mykit.core.constants import EV2RY
from mykit.core.dos import Dos, DosError

goodEgrid = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
badEgrid = [-4, -3, -2, -1, 0, 1]
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
badDos1 = [1.0, 2.0, 3.0, 4.0],
badDos2 = [
    [3.0, 4.0],
    [1.0, 2.0],
    [0.0, 0.0],
    [0.0, ]
]
badDos = (badDos1, badDos2)
nedos, nspins = np.shape(goodDos)


class test_dos_initialize(ut.TestCase):

    def test_raise_for_inconsistent_egrid_dos(self):
        for b in badDos:
            self.assertRaisesRegex(DosError, r"Inconsistent shape: *",
                                   Dos, goodEgrid, b, 0.0)
        self.assertRaisesRegex(DosError, r"Inconsistent shape: *",
                               Dos, badEgrid, goodDos, 0.0)

    def test_properties(self):
        dos = Dos(goodEgrid, goodDos, efermi=1.0, unit="ev")
        self.assertFalse(dos.hasProjection)
        self.assertEqual(dos.nedos, nedos)
        self.assertEqual(dos.nspins, nspins)
        self.assertEqual('ev', dos.unit)
        self.assertTrue(np.allclose(goodEgrid, dos.edos))
        self.assertTrue(np.allclose(goodDos, dos.dos))
        dos.unit = "ry"
        self.assertAlmostEqual(dos.efermi, EV2RY)
        # None for atoms, projs and pDos when no projection was parsed
        self.assertEqual(None, dos.atoms)
        self.assertEqual(0, dos.natoms)
        self.assertEqual(None, dos.projs)
        self.assertEqual(0, dos.nprojs)
        self.assertEqual(None, dos.pDos)

    def test_sum_proj(self):
        dos = Dos(goodEgrid, goodDos, efermi=0.0, unit="ev")
        # no projection was parsed. empty list
        self.assertListEqual([], dos._get_atom_indices(None))
        self.assertListEqual([], dos._get_proj_indices(None))
        self.assertTrue(np.array_equal(dos.sum_atom_proj_comp(
            fail_one=True), np.ones((nedos, nspins))))
        # raise for bad fail_one type
        self.assertRaisesRegex(TypeError, r"fail_one should be bool type.",
                               dos.sum_atom_proj_comp, fail_one=3)


if __name__ == '__main__':
    ut.main()
