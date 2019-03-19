#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.core.utils import get_matched_files
from mykit.vasp.xml import Vasprunxml, VasprunxmlError


class test_vasprunxml_read(ut.TestCase):

    def test_scf_xml(self):
        '''Test reading XMLs for SCF calculations (LORBIT not set)
        '''
        dataDir = 'vasprun_scf'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            vxml = Vasprunxml(fn)
            typeMapping = vxml.typeMapping
            # get all index
            self.assertListEqual(list(range(vxml.natoms)), \
                vxml.get_atom_index())
            self.assertFalse(vxml.dosGrid is None)
            self.assertFalse(vxml.totalDos is None)
            self.assertFalse(vxml.dos is None)
            vxml.get_atom_index(0)
            vxml.get_atom_index(-1)
            vxml.get_atom_index(typeMapping[0])
            self.assertRaisesRegex(VasprunxmlError, \
                r"Atom type not found: *", \
                vxml.get_atom_index, "UNKNOWNSYMBOL")
            # empty properties
            self.assertTrue(vxml.projs is None)
            self.assertEqual(0, vxml.nprojs)
            self.assertTrue(vxml.pDos is None)
            self.assertTrue(vxml.pWave is None)

    def test_band_xml(self):
        dataDir = 'vasprun_band'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            _vxml = Vasprunxml(fn)

    def test_opt_xml(self):
        dataDir = 'vasprun_opt'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            _vxml = Vasprunxml(fn)

    def test_pdos_xml(self):
        dataDir = 'vasprun_partial'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            vxml = Vasprunxml(fn)
            self.assertFalse(vxml.pDos is None)
            bs = vxml.load_band()
            bs.pWave
            bs.atoms
            bs.projs
            self.assertTrue(bs.hasProjection)
                

if __name__ == '__main__':
    ut.main()
