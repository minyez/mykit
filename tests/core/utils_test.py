#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut

from mykit.core.utils import get_dirpath, trim_after, trim_before


class file_and_path(ut.TestCase):

    def test_get_dirpath(self):
        self.assertEqual("/usr/bin", get_dirpath("/usr/bin"))
        self.assertEqual("/bin", get_dirpath("/bin/ls"))


class test_string_manipulation(ut.TestCase):

    def test_trim_after(self):
        self.assertEqual("abc", trim_after("abc#defg", r'#'))
        self.assertEqual("I have Fortran", trim_after("I have Fortran!comment", r'!'))
        self.assertEqual("Fe", trim_after("Fe1", r'\d'))
        self.assertEqual("Cd", trim_after("Cd2Fe", r'\d'))
    
    def test_trim_before(self):
        self.assertEqual("#defg", trim_before("abc#defg", r'#'))
        self.assertEqual("!comment", trim_before("I have Fortran!comment", r'!'))
        self.assertEqual("1", trim_before("Fe1", r'\d'))
        self.assertEqual("2Fe", trim_before("Cd2Fe", r'\d'))


if __name__ == "__main__":
    ut.main()
