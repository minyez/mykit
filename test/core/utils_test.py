#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut

import numpy as np

from mykit.core.utils import (Cif, conv_estimate_number, conv_string,
                              find_vol_dirs, get_dirpath, get_file_ext,
                              get_latt_vecs_from_latt_consts, get_str_indices,
                              get_str_indices_by_iden, trim_after, trim_before,
                              trim_both_sides)


class file_and_path(ut.TestCase):

    def test_get_dirpath(self):
        self.assertEqual("/usr/bin", get_dirpath("/usr/bin"))
        self.assertEqual("/bin", get_dirpath("/bin/ls"))

    def test_find_vol_dirs(self):
        # no volume directories at cwd
        self.assertListEqual([], find_vol_dirs())

    def test_get_file_ext(self):
        self.assertEqual("cif", get_file_ext('abc.cif'))
        self.assertEqual("struct", get_file_ext('a/b/c/c.struct'))
        self.assertEqual("", get_file_ext('a/b/c/c'))
        self.assertEqual(None, get_file_ext(get_dirpath(__file__)))


class test_string_manipulation(ut.TestCase):

    def test_trim_after(self):
        self.assertEqual("abc", trim_after("abc#defg", r'#'))
        self.assertEqual("abc#", trim_after("abc#defg", r'#', include_pattern=True))
        self.assertEqual("I have Fortran", \
            trim_after("I have Fortran!comment", r'!'))
        self.assertEqual("I have Fortran!", \
            trim_after("I have Fortran!comment", r'!', include_pattern=True))
        self.assertEqual("Fe", trim_after("Fe1", r'\d'))
        self.assertEqual("Cd2", trim_after("Cd2Fe", r'\d', include_pattern=True))
    
    def test_trim_before(self):
        self.assertEqual("defg", trim_before("abc#defg", r'#'))
        self.assertEqual("#defg", trim_before("abc#defg", r'#', include_pattern=True))
        self.assertEqual("comment", trim_before("I have Fortran!comment", r'!'))
        self.assertEqual("!comment", \
            trim_before("I have Fortran!comment", r'!', include_pattern=True))
        self.assertEqual("P", trim_before("Fe1P", r'\d'))
        self.assertEqual("1P", trim_before("Fe1P", r'\d', include_pattern=True))
        self.assertEqual("Fe", trim_before("Cd2Fe", r'\d'))
        self.assertEqual("2Fe", trim_before("Cd2Fe", r'\d', include_pattern=True))

    def test_trim_both(self):
        string = "WFFIL  EF=0.9725 (WFFIL, WFPRI, ENFIL, SUPWF)"
        self.assertEqual("0.9725 ", \
            trim_both_sides(string, r"=", r"\("))
        self.assertEqual("=0.9725 (", \
            trim_both_sides(string, r"=", r"\(", include_pattern=True))

    def test_conv_string(self):
        string = "1; ABC=5. "
        # unsupported conversion type
        self.assertRaises(AssertionError, conv_string, string, list, 0)
        # single value
        self.assertEqual(1, \
            conv_string(string, int, 0, strips=";"))
        # multiple values
        self.assertListEqual([1, 5], \
            conv_string(string, int, 0, 2, sep=r"[;=]", strips="."))
        # # inpropriate value to convert
        # self.assertRaises(ValueError, conv_string, \
        #     string, int, 0, 2, sep=r"[;=]")

    def test_conv_estimate_number(self):
        self.assertEqual(3.23, conv_estimate_number('3.2(3)'))


class test_indicies_searching(ut.TestCase):

    container = ['A', 'b', 'C', 'M', 'A', 'D', 'b']

    def test_get_str_indices(self):
        self.assertListEqual([0, 4], get_str_indices(self.container, 'A'))
        self.assertListEqual([5,], get_str_indices(self.container, 'D'))
        self.assertListEqual([], get_str_indices(self.container, 'Xz'))

    def test_get_str_indices_by_iden(self):
        self.assertListEqual([0, 4, 1], \
            get_str_indices_by_iden(self.container, ['A', 1]))
        self.assertListEqual([1, 6, 0], \
            get_str_indices_by_iden(self.container, ['b', 0]))
        self.assertListEqual([1, 6], \
            get_str_indices_by_iden(self.container, ['b', 20]))


class test_structure_related(ut.TestCase):
    '''Test utilities related with structure and cell creation
    '''

    def test_decode_equiv_pos_string(self):
        self.assertRaisesRegex(ValueError, \
            r"s does not seem to be a symmetry operation string", 
            Cif.decode_equiv_pos_string, 'x, y, z, a')
        self.assertRaisesRegex(ValueError, \
            r"s does not seem to be a symmetry operation string", 
            Cif.decode_equiv_pos_string, 'x, y,')
    
    def test_get_latt_vecs_from_latt_consts(self):
        a = 2.0
        b = 3.0
        c = 4.0
        # alpha = 60
        # beta = 60
        # gamma = 60
        self.assertTrue(np.allclose(\
            [[a,0,0],[0,b,0],[0,0,c]], \
                get_latt_vecs_from_latt_consts(a, b, c)))


if __name__ == "__main__":
    ut.main()
