#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
import numpy as np
from mykit.core.constants import au2ang, ang2au
from mykit.core.lattice import lattice, latticeError

class simple_cubic_lattice(unittest.TestCase):

    _a = 5.0
    _cell = [[_a, 0.0, 0.0],
             [0.0, _a, 0.0],
             [0.0, 0.0, _a]]
    _atoms = ["C"]
    _frac = 0.5
    _pos = [[_frac, 0.0, 0.0]]

    _latt = lattice(_cell, _atoms, _pos, unit="ang", coordSys="D")
    
    def test_read_only_properties(self):
        self.assertAlmostEqual(pow(self._a, 3), self._latt.vol)
        try:
            self._latt.vol = 10.0
        except AttributeError as _err:
            self.assertIn("can't set attribute", _err.args)
        else:
            self.assertTrue(False, "lattice.vol should be read-only")
        self.assertEqual(1, self._latt.natoms)
        try:
            self._latt.natoms = 2
        except AttributeError as _err:
            self.assertIn("can't set attribute", _err.args)
        else:
            self.assertTrue(False, "lattice.natoms should be read-only")
        self.assertTupleEqual(tuple(self._atoms), self._latt.atomTypes)
        try:
            self._latt.atomTypes = ["B"]
        except AttributeError as _err:
            self.assertIn("can't set attribute", _err.args)
        else:
            self.assertTrue(False, "lattice.atomTypes should be read-only")
        self.assertTupleEqual(tuple(range(len(self._atoms))), self._latt.typeIndex)
        try:
            self._latt.typeIndex = (0)
        except AttributeError as _err:
            self.assertIn("can't set attribute", _err.args)
        else:
            self.assertTrue(False, "lattice.typeIndex should be read-only")

    def test_magic(self):
        self.assertEqual(1, len(self._latt))
        self.assertTupleEqual(tuple(self._pos[0]), tuple(self._latt[0]))

    def test_unit_conv(self):
        # ang2au
        self._latt.unit = 'au'
        _cell = self._latt.get_latt()[0]
        if self._latt._dtype == 'float32':
            self.assertAlmostEqual(_cell[0,0], self._a * ang2au, places=5)
        else:
            self.assertAlmostEqual(_cell[0,0], self._a * ang2au)
        # au2ang
        self._latt.unit = 'ang'
        _cell = self._latt.get_latt()[0]
        self.assertEqual(_cell[0,0], self._a)

    def test_coord_conv(self):
        # direct2cart
        self._latt.coordSys = 'C'
        self.assertTupleEqual(tuple(self._latt[0]), (self._frac * self._a, 0.0, 0.0))
        # cart2direct
        self._latt.coordSys = 'D'
        self.assertEqual(self._latt[0][0], self._frac)
        

class lattice_raise(unittest.TestCase):
    
    def test_bad_cell(self):
        _cell = [[5.0, 0.0, 0.0],
                 [0.0, 0.0, 5.0]]
        _atoms = ["C"]
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(latticeError, lattice, _cell, _atoms, _pos)
        _cell = [[5.0, 0.0],
                 [0.0, 5.0, 0.0],
                 [0.0, 0.0, 5.0]]
        self.assertRaises(latticeError, lattice, _cell, _atoms, _pos)

    def test_bad_atoms_pos(self):
        _cell = [[5.0, 0.0, 0.0],
                 [0.0, 5.0, 0.0],
                 [0.0, 0.0, 5.0]]

        _atoms = []
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(latticeError, lattice, _cell, _atoms, _pos)

        _atoms = ["C"]
        _pos = [0.0, 0.0, 0.0]
        self.assertRaises(latticeError, lattice, _cell, _atoms, _pos)

        _atoms = ["C"]
        _pos = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]
        self.assertRaises(latticeError, lattice, _cell, _atoms, _pos)

        _atoms = ["C", "C"]
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(latticeError, lattice, _cell, _atoms, _pos)


class lattice_factory_method(unittest.TestCase):
    '''Test the class methods to generate commonly used lattice structure'''

    def test_bravis_cubic(self):
        _pc = lattice.bravis_cP("C", aLatt=5.0, coordSys="D")
        self.assertEqual(1, len(_pc))
        self.assertEqual("D", _pc.coordSys)
        _bcc = lattice.bravis_cI("C", aLatt=5.0, unit="au")
        self.assertEqual("au", _bcc.unit)
        self.assertEqual(2, len(_bcc))
        _fcc = lattice.bravis_cF("C", aLatt=5.0)
        self.assertEqual(4, len(_fcc))


class lattice_element_utils(unittest.TestCase):
    '''Test the utils to manipulate element lists for ``lattice`` use
    '''
    pass

if __name__ == "__main__":
    unittest.main()
    
