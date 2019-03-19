#!/usr/bin/env python3
# coding = utf-8

import json
import os
import re
import unittest as ut

import numpy as np

from mykit.core.bandstructure import BandStructure as BS
from mykit.core.bandstructure import BandStructureError as BSE
from mykit.core.bandstructure import _check_eigen_occ_consistency

goodEigen = [
    [
        [1,2,3],
        [11,22,33],
    ],
]
badEigen = [[100,102,103], [210,212,213]]
goodOcc = [
    [
        [1.0,1.0,0.0],
        [1.0,1.0,0.0],
    ],
]
badOcc = [[1.0,1.0,0.0],[1.0,1.0,0.0]]
nspins, nkpts, nbands = np.shape(goodEigen)

class test_check_consistency(ut.TestCase):

    def test_check_eigen_occ(self):
        self.assertTupleEqual(_check_eigen_occ_consistency(goodEigen, goodOcc), \
            (nspins, nkpts, nbands))
        self.assertTupleEqual(_check_eigen_occ_consistency(badEigen, badOcc), \
            ())


class test_BS_no_projection(ut.TestCase):

    def test_raise_inconsistent_eigen_occ(self):
        self.assertRaisesRegex(BSE, r"Bad eigen and occ shapes *", \
            BS, badEigen, goodOcc)
        self.assertRaisesRegex(BSE, r"Bad eigen and occ shapes *", \
            BS, goodEigen, badOcc)

    def test_properties(self):
        bs = BS(goodEigen, goodOcc)
        self.assertTrue(np.allclose(goodEigen, bs.eigen))
        self.assertTrue(np.allclose(goodOcc, bs.occ))
        self.assertEqual(bs.nspins, nspins)
        self.assertEqual(bs.nkpts, nkpts)
        self.assertEqual(bs.nbands, nbands)
        # empty properties when initialized without projections
        self.assertEqual(None, bs.atoms)
        self.assertEqual(None, bs.projs)
        self.assertEqual(None, bs.pWave)
        self.assertFalse(bs.hasProjection)

class test_BS_projection(ut.TestCase):

    def test_reading_in_good_projection(self):
        dataDirPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
            '..', 'testdata')
        countGood = 0
        for fn in os.listdir(dataDirPath):
            if re.match(r'Bandstructure_proj_[\d]+\.json', fn):
                path = os.path.join(dataDirPath, fn)
                with open(path, 'r') as f:
                    j = json.load(f)
                shape = j.pop("shape")
                # self.assertTupleEqual(np.shape(j["pWave"]), tuple(j["shape"]))
                self.assertTupleEqual(np.shape(j["pWave"]), tuple(shape))
                eigen = np.random.random(shape[:3])
                occ = np.random.choice([0.0, 1.0], shape[:3])
                bs = BS(eigen, occ, projected=j)
                self.assertListEqual(bs.atoms, j["atoms"])
                self.assertListEqual(bs.projs, j["projs"])
                self.assertTrue(bs.hasProjection)
                countGood += 1
        print("Processed {} good band structure projections".format(countGood))


if __name__ == '__main__':
    ut.main()
