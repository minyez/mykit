#!/usr/bin/env python3
# coding = utf-8

import json
import os
import re
import unittest as ut

import numpy as np

from mykit.core.bandstructure import BandStructure as BS
from mykit.core.bandstructure import BandStructureError as BSE
from mykit.core.bandstructure import _check_eigen_occ_weight_consistency
from mykit.core.utils import get_matched_files

goodEigen = [
    [
        [1,2,3],
        [11,22,33],
    ],
]
badEigen = [[100,102,103], [210,212,213]]
goodOcc = [
    [
        [2.0,2.0,0.0],
        [2.0,2.0,0.0],
    ],
]
badOcc = [[1.0,1.0,0.0],[1.0,1.0,0.0]]
goodWeight = [1, 4]
badWeight = [1,]
nspins, nkpts, nbands = np.shape(goodEigen)

class test_check_consistency(ut.TestCase):

    def test_check_eigen_occ(self):
        self.assertTupleEqual(_check_eigen_occ_weight_consistency(goodEigen, goodOcc, goodWeight), \
            (nspins, nkpts, nbands))
        self.assertTupleEqual(_check_eigen_occ_weight_consistency(badEigen, badOcc, goodWeight), \
            ())
        self.assertTupleEqual(_check_eigen_occ_weight_consistency(goodEigen, goodOcc, badWeight), \
            ())


class test_BS_no_projection(ut.TestCase):

    def test_raise_inconsistent_eigen_occ(self):
        self.assertRaisesRegex(BSE, r"Bad eigen, occ and weight shapes *", \
            BS, badEigen, goodOcc, goodWeight)
        self.assertRaisesRegex(BSE, r"Bad eigen, occ and weight shapes *", \
            BS, goodEigen, badOcc, goodWeight)
        self.assertRaisesRegex(BSE, r"Bad eigen, occ and weight shapes *", \
            BS, goodEigen, goodOcc, badWeight)

    def test_properties(self):
        bs = BS(goodEigen, goodOcc, goodWeight)
        self.assertTrue(np.allclose(goodEigen, bs.eigen))
        self.assertTrue(np.allclose(goodOcc, bs.occ))
        self.assertEqual(bs.nspins, nspins)
        self.assertEqual(bs.nkpts, nkpts)
        self.assertEqual(bs.nbands, nbands)
        self.assertTupleEqual((nspins, nkpts), np.shape(bs.ivbm))
        self.assertTupleEqual((nspins, nkpts), np.shape(bs.icbm))
        self.assertTupleEqual((nspins, nkpts), np.shape(bs.vbm))
        self.assertTupleEqual((nspins, nkpts), np.shape(bs.cbm))
        self.assertTupleEqual((nspins, nkpts), np.shape(bs.directGap))
        self.assertTupleEqual((nspins,), np.shape(bs.fundGap))
        self.assertTupleEqual((nspins, 2), np.shape(bs.fundTrans))
        self.assertTupleEqual((nspins,), np.shape(bs.kAvgGap))
        # empty properties when initialized without projections
        self.assertTrue(bs.atoms is None)
        self.assertTrue(bs.projs is None)
        self.assertTrue(bs.pWave is None)
        self.assertFalse(bs.hasProjection)


class test_BS_projection(ut.TestCase):

    dataDirPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
        '..', 'testdata')

    def test_reading_in_good_projection(self):
        countGood = 0
        for path in get_matched_files(self.dataDirPath, r'Bandstructure_proj_[\d]+\.json'):
            with open(path, 'r') as f:
                j = json.load(f)
            shape = j.pop("shape")
            # self.assertTupleEqual(np.shape(j["pWave"]), tuple(j["shape"]))
            self.assertTupleEqual(np.shape(j["pWave"]), tuple(shape))
            # as occ is randomized, there are cases that warnings infinity CBM
            # it is totally okay.
            eigen = np.random.random(shape[:3])
            occ = np.random.choice([0.0, 1.0], shape[:3])
            weight = np.random.random(shape[1])
            bs = BS(eigen, occ, weight, projected=j)
            self.assertListEqual(bs.atoms, j["atoms"])
            self.assertListEqual(bs.projs, j["projs"])
            self.assertTrue(bs.hasProjection)
            countGood += 1
        print("Processed {} good band structure projections".format(countGood))




if __name__ == '__main__':
    ut.main()
