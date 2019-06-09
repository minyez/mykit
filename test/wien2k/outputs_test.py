#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.wien2k.outputs import Vsp

class test_vsp(ut.TestCase):

    def read(self):
        dataDir = 'vsp'
        dataDirPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
            '..', 'testdata', 'wien2k', dataDir)
        for fn in os.listdir(dataDirPath):
            path = os.path.join(dataDirPath, fn)
            _v = Vsp(path)


if __name__ == "__main__":
    ut.main()
