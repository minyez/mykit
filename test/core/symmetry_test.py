#!/usr/bin/env python3
# coding = utf-8

import unittest as ut

from mykit.core.metadata._spk import (_special_kpoints, cond_a_lt_b,
                                      cond_abc_invsq, cond_any, cond_c_lt_a,
                                      cond_curt_a_lt_sqrt_c,
                                      cond_max_bc_noneq_left)
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

class test_spk_dict(ut.TestCase):
    '''Test the integrity of _special_kpoints dictionary
    '''
    def test_cond_func(self):
        self.assertEqual(cond_any(1,2,3), 0)

        self.assertEqual(cond_a_lt_b(2,1,3), 0)
        self.assertEqual(cond_a_lt_b(1,2,3), 1)
        self.assertRaisesRegex(ValueError, r"Fail to determine special kpoints set", cond_a_lt_b, 1, 1, 3)

        self.assertEqual(cond_c_lt_a(2,1,3), 0)
        self.assertEqual(cond_c_lt_a(3,1,2), 1)
        self.assertRaisesRegex(ValueError, r"Fail to determine special kpoints set", cond_c_lt_a, 1, 3, 1)

        self.assertEqual(cond_curt_a_lt_sqrt_c(8,1,3), 0)
        self.assertEqual(cond_curt_a_lt_sqrt_c(8,1,6), 1)
        self.assertRaisesRegex(ValueError, r"Fail to determine special kpoints set", cond_curt_a_lt_sqrt_c, 8, 1, 4)

        self.assertEqual(cond_max_bc_noneq_left(2,5,3), 0)
        self.assertEqual(cond_max_bc_noneq_left(3,5,2), 0)
        self.assertEqual(cond_max_bc_noneq_left(2,3,5), 1)
        self.assertEqual(cond_max_bc_noneq_left(3,2,5), 1)
        self.assertRaisesRegex(ValueError, r"Fail to determine special kpoints set", cond_max_bc_noneq_left, 2, 5, 2)
        self.assertRaisesRegex(ValueError, r"Fail to determine special kpoints set", cond_max_bc_noneq_left, 2, 2, 5)
        self.assertRaisesRegex(ValueError, r"Fail to determine special kpoints set", cond_max_bc_noneq_left, 5, 2, 2)

        self.assertEqual(cond_abc_invsq(0.1,0.5,0.5), 0)
        self.assertEqual(cond_abc_invsq(0.5,0.5,0.1), 1)
        self.assertEqual(cond_abc_invsq(0.2,0.2,0.5), 2)
        self.assertRaisesRegex(ValueError, r"Fail to determine special kpoints set", cond_abc_invsq, 0.5, 0.1, 0.5)


    def test_number_of_spg(self):
        for i in range(1, 231):
            self.assertIn(i, _special_kpoints.keys())
        self.assertEqual(230, len(_special_kpoints.keys()))

    def test_has_cond_and_spPrim(self):
        for k, v in _special_kpoints.items():
            self.assertIn("cond", v, msg="bad key {}".format(k))
            self.assertIn("spPrim", v, msg="bad key {}".format(k))
            self.assertTrue(isinstance(v["spPrim"], list), msg="bad key {}".format(k))

    def test_number_of_sets(self):
        for k, v in _special_kpoints.items():
            if v["cond"] == cond_any:
                self.assertEqual(1, len(v["spPrim"]), msg="bad key {}".format(k))
            if v["cond"] == cond_abc_invsq:
                self.assertEqual(3, len(v["spPrim"]), msg="bad key {}".format(k))
            hasTwoSets = [
                cond_a_lt_b, 
                cond_c_lt_a, 
                cond_curt_a_lt_sqrt_c,
                cond_max_bc_noneq_left,
                ]
            if v["cond"] in hasTwoSets:
                self.assertEqual(2, len(v["spPrim"]), msg="bad key {}".format(k))


if __name__ == '__main__':
    ut.main()
