#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut
from mykit.core.utils import get_dirpath, trim_comment

class file_and_path(ut.TestCase):

    def test_get_dirpath(self):
        self.assertEqual("/usr/bin", get_dirpath("/usr/bin"))
        self.assertEqual("/bin", get_dirpath("/bin/ls"))


class test_string_manipulation(ut.TestCase):

    def test_trim_comment(self):
        self.assertEqual("abc", trim_comment("abc#defg", r'#'))
        self.assertEqual("I have Fortran", trim_comment("I have Fortran!comment", r'!'))

if __name__ == "__main__":
    ut.main()