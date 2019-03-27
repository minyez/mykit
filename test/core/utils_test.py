#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut

from mykit.core.utils import (conv_string, find_vol_dirs, get_dirpath,
                              get_file_ext, get_str_indices,
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
        self.assertEqual("struct", get_file_ext('a/b/c/c'))
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


if __name__ == "__main__":
    ut.main()
