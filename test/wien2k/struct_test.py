#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.wien2k.struct import Struct

class test_struct(ut.TestCase):

    def test_read_from_file(self):
        dataDir = 'struct'
        dataDirPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
            '..', 'testdata', 'wien2k', dataDir)
        for fn in os.listdir(dataDirPath):
            path = os.path.join(dataDirPath, fn)
            _s = Struct.read_from_file(path)


if __name__ == '__main__':
    ut.main()