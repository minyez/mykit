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
    
    def test_properties(self):
        self.assertAlmostEqual(pow(self._a, 3), self._latt.vol)
        try:
            self._latt.vol = 10.0
        except AttributeError as _err:
            self.assertIn("can't set attribute", _err.args)
        else:
            self.assertTrue(False, "lattice.vol is not read-only")
        self.assertEqual(1, self._latt.natoms)
        try:
            self._latt.natoms = 2
        except AttributeError as _err:
            self.assertIn("can't set attribute", _err.args)
        else:
            self.assertTrue(False, "lattice.natoms is not read-only")

    def test_magic(self):
        self.assertEqual(1, len(self._latt))
        self.assertTupleEqual(tuple(self._pos[0]), tuple(self._latt[0]))

    def test_unit_conv(self):
        # ang2au
        self._latt.unit = 'au'
        if self._latt._dtype == 'float32':
            self.assertAlmostEqual(self._latt.cell[0,0], self._a * ang2au, places=5)
        else:
            self.assertAlmostEqual(self._latt.cell[0,0], self._a * ang2au)
        # au2ang
        self._latt.unit = 'ang'
        self.assertEqual(self._latt.cell[0,0], self._a)

    def test_coord_conv(self):
        # direct2cart
        self._latt.coordSys = 'C'
        self.assertTupleEqual(tuple(self._latt.pos[0]), (self._frac * self._a, 0.0, 0.0))
        # cart2direct
        self._latt.coordSys = 'D'
        self.assertEqual(self._latt.pos[0,0], self._frac)
        

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


if __name__ == "__main__":
    unittest.main()
    
