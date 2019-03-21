#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

import numpy as np

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
            vxml.ntypes
            vxml.natomsPerType
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
        '''Test reading XMLs for band calculations (LORBIT set or not)
        '''
        dataDir = 'vasprun_band'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            vxml = Vasprunxml(fn)
            self.assertEqual(vxml.kmode, "L")
            self.assertTupleEqual(np.shape(vxml.weight), (vxml.nibzkpt,))
            self.assertTupleEqual(np.shape(vxml.kpoints), (vxml.nibzkpt, 3))
            self.assertTupleEqual(np.shape(vxml.kptsWeight), (vxml.nibzkpt, 4))
            bs = vxml.load_band()
            self.assertAlmostEqual(bs.nelect, vxml.nelect, places=4)
            self.assertTrue(bs.hasKvec)
            self.assertTrue(bs.isKpath)
            bs.kvec
            bsTrimed = vxml.load_band(1)
            self.assertEqual(1, bs.nkpts - bsTrimed.nkpts)

    def test_mixed_k_band_xml(self):
        '''Test reading XMLs for band calculations with manual input kpoints
        in case of SCAN and HF band calculations
        '''
        dataDir = 'vasprun_mixed_k_band'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            vxml = Vasprunxml(fn)
            bsMix = vxml.load_band()
            bsBand = vxml.load_band(kTrimBefore=20)
            self.assertEqual(bsMix.nkpts - bsBand.nkpts, 20)
            self.assertTrue(np.allclose(bsBand.weight, np.ones(bsBand.nkpts)))
            self.assertTrue(bsBand.isKpath)


    def test_opt_xml(self):
        '''Test reading XMLs for geometry optimization
        '''
        dataDir = 'vasprun_opt'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            _vxml = Vasprunxml(fn)

    def test_pdos_xml(self):
        '''Test reading XMLs with LORBIT set
        '''
        dataDir = 'vasprun_partial'
        dataDirPath = os.path.join(os.path.dirname(__file__), '..', \
            'testdata', 'vasp', dataDir)
        for fn in get_matched_files(dataDirPath, r"vasprun*"):
            vxml = Vasprunxml(fn)
            self.assertFalse(vxml.pDos is None)
            bs = vxml.load_band()
            self.assertAlmostEqual(bs.nelect, vxml.nelect, places=4)
            self.assertTrue(bs.hasProjection)
            # bs.pWave
            # bs.atoms
            # bs.projs
                

if __name__ == '__main__':
    ut.main()
