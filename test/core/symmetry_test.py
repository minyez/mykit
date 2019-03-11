#!/usr/bin/env python3
# coding = utf-8

import unittest as ut

from mykit.core.symmetry import (Symmetry, SymmetryError, space_group,
                                 special_kpoints)


class test_space_group(ut.TestCase):

    def test_get_spg_symbol(self):
        self.assertEqual(space_group.get_spg_symbol(1), 'P1')
        self.assertEqual(space_group.get_spg_symbol(22), 'F222')
        self.assertEqual(space_group.get_spg_symbol(64), 'Cmce')
        self.assertEqual(space_group.get_spg_symbol(123), 'P4mmm')
        self.assertEqual(space_group.get_spg_symbol(161), 'R3c')
        self.assertEqual(space_group.get_spg_symbol(208), 'P4232')
        self.assertEqual(space_group.get_spg_symbol(230), 'Ia-3d')
        # test raises
        self.assertRaisesRegex(SymmetryError, r"string received. id should be int", \
            space_group.get_spg_symbol, 'P1')
        self.assertRaisesRegex(SymmetryError, r"Invalid space group id \(1~230\): 231", \
            space_group.get_spg_symbol, 231)
        self.assertRaisesRegex(SymmetryError, r"Invalid space group id \(1~230\): -1", \
            space_group.get_spg_symbol, -1)
        


if __name__ == '__main__':
    ut.main()
