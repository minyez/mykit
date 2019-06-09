#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.wien2k.utils import get_casename, find_complex_file, get_default_r0, get_default_rmt, get_z

class test_get_elements_info(ut.TestCase):

    def test_r0(self):
        self.assertEqual(get_default_r0('Ne1'), get_default_r0('Ne10'))

    def test_rmt(self):
        self.assertEqual(get_default_rmt('Mg2'), get_default_rmt('Mg'))
    
    def test_z(self):
        self.assertEqual(get_z('X'), 0.0)
        self.assertEqual(get_z('Pb'), 82.0)
        self.assertEqual(get_z('Cl'), 17.0)
        self.assertEqual(get_z('Cu'), 29.0)
        self.assertEqual(get_z('Ag'), 47.0)


class test_search_file(ut.TestCase):

    def test_find_complex_file(self):

        self.assertRaises(FileNotFoundError, find_complex_file, 'fake_case', 'in')
        self.assertEqual(get_casename(os.path.dirname(__file__)), 'wien2k')


if __name__ == '__main__':
    ut.main()