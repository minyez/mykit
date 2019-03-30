#!/usr/bin/env python3
# coding=utf-8

import unittest as ut
from collections import OrderedDict

from mykit.core.elements import ATOMIC_WEIGHT, ELEMENT_SYMBOLS


class test_elem_symbols(ut.TestCase):

    def test_consistency(self):
        self.assertEqual(len(ATOMIC_WEIGHT), len(ELEMENT_SYMBOLS))
        # check duplicates in element symbols
        od = OrderedDict.fromkeys(ELEMENT_SYMBOLS)
        self.assertEqual(len(ELEMENT_SYMBOLS), len(od.keys()))

    def test_right_order(self):
        oi = -1
        for e in ['H', 'C', 'Na', 'S', 'Ga', 'As', 'Pd', 'La']:
            self.assertIn(e, ELEMENT_SYMBOLS)
            ni = ELEMENT_SYMBOLS.index(e)
            self.assertTrue(ni > oi)
            oi = ni


if __name__ == '__main__':
    ut.main()
