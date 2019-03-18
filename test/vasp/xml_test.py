#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.vasp.xml import Vasprunxml, VasprunxmlError


class test_vasprunxml_read(ut.TestCase):

    def test_static_xml(self):
        dataDir = 'vasprun_static'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', dataDir)
        read_all_xml(dataDirPath)

    def test_opt_xml(self):
        dataDir = 'vasprun_opt'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', dataDir)
        read_all_xml(dataDirPath)


def read_all_xml(dirpath):        
    for d in os.listdir(dirpath):
        path = os.path.join(dirpath, d)
        _vxml = Vasprunxml(path)


if __name__ == '__main__':
    ut.main()

