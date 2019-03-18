#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut

from mykit.core.utils import (get_dirpath, trim_after, trim_before,
                              trim_both_sides)


class file_and_path(ut.TestCase):

    def test_get_dirpath(self):
        self.assertEqual("/usr/bin", get_dirpath("/usr/bin"))
        self.assertEqual("/bin", get_dirpath("/bin/ls"))


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


if __name__ == "__main__":
    ut.main()
