#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut
from fnmatch import fnmatch

from mykit.vasp.xml import Vasprunxml, VasprunxmlError


class test_vasprunxml_read(ut.TestCase):

    def test_static_xml(self):
        dataDir = 'vasprun_static'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in os.listdir(dataDirPath):
            if fnmatch(fn, r"vasprun*"):
                path = os.path.join(dataDirPath, fn)
                _vxml = Vasprunxml(path)
                _vxml.get_atom_index()
                _vxml.get_atom_index(0)
                _vxml.get_atom_index(-1)

    def test_opt_xml(self):
        dataDir = 'vasprun_opt'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in os.listdir(dataDirPath):
            if fnmatch(fn, r"vasprun*"):
                path = os.path.join(dataDirPath, fn)
                _vxml = Vasprunxml(path)
                _vxml.get_atom_index()
                _vxml.get_atom_index(0)
                _vxml.get_atom_index(-1)

    def test_pdos_xml(self):
        dataDir = 'vasprun_partial'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in os.listdir(dataDirPath):
            if fnmatch(fn, r"vasprun*"):
                path = os.path.join(dataDirPath, fn)
                _vxml = Vasprunxml(path)
                _vxml.get_atom_index()
                _vxml.get_atom_index(0)
                _vxml.get_atom_index(-1)
                _vxml.load_band()


if __name__ == '__main__':
    ut.main()
