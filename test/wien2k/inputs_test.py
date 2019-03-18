#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.wien2k.inputs import In1, InputError


class test_In1(ut.TestCase):

    def test_read_from_file(self):
        dataDir = 'in1'
        dataDirPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
            '..', 'testdata', 'wien2k', dataDir)
        for fn in os.listdir(dataDirPath):
            path = os.path.join(dataDirPath, fn)
            In1.read_from_file(path)

if __name__ == '__main__':
    ut.main()
