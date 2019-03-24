#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

import numpy as np

from mykit.core.utils import get_matched_files
from mykit.vasp.incar import Incar
from mykit.vasp.outcar import Outcar, OutcarError, get_value


class test_Outcar(ut.TestCase):

    dataDir = os.path.join(os.path.dirname(__file__), '..', 'testdata', 'vasp')

    def test_reading_and_properties(self):

        for fn in get_matched_files(self.dataDir, r"OUTCAR_[0-9]+"):
            oc = Outcar(fn)
            # check outcar reading and initialiaztion
            self.assertEqual(len(oc.kpoints), oc.nkpts)
            self.assertEqual(len(oc.weight), oc.nkpts)
            self.assertEqual(len(oc._iterations), oc.nIonSteps)
            self.assertIsInstance(oc.incar, Incar)
            # test load band structure
            bs = oc.load_band()
            self.assertTupleEqual((oc.nspins, oc.nkpts, oc.nbands), \
                np.shape(bs.eigen))
            self.assertEqual(bs.nbands, oc.nbands)
            self.assertEqual(bs.efermi, oc.efermi)

    def test_utilities(self):
        for fn in get_matched_files(self.dataDir, r"OUTCAR_[0-9]+"):
            _fermi = get_value("fermi", fn)
            _encut = get_value("ENCUT", fn)
            



if __name__ == '__main__':
    ut.main()
