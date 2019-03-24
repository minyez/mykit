#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.core.utils import get_matched_files
from mykit.vasp.outcar import Outcar, OutcarError


class test_Outcar(ut.TestCase):

    dataDir = os.path.join(os.path.dirname(__file__), '..', 'testdata', 'vasp')

    def test_reading(self):
        # read good OUTCARs
        for fn in get_matched_files(self.dataDir, r"OUTCAR_[0-9]+"):
            oc = Outcar(fn)
            self.assertEqual(len(oc.kpoints), oc.nkpts)
            self.assertEqual(len(oc.weight), oc.nkpts)


if __name__ == '__main__':
    ut.main()
