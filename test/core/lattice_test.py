#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut

import numpy as np

from mykit.core.constants import ang2au, au2ang
from mykit.core.lattice import (LatticeError, atoms_from_sym_nat, axis_list,
                                lattice, periodic_duplicates_in_cell,
                                sym_nat_from_atoms)


class simple_cubic_lattice(ut.TestCase):

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
        self.assertListEqual(list(self._atoms), self._latt.atomTypes)
        try:
            self._latt.atomTypes = ["B"]
        except AttributeError as _err:
            self.assertIn("can't set attribute", _err.args)
        else:
            self.assertTrue(False, "lattice.atomTypes should be read-only")
        self.assertListEqual(list(range(len(self._atoms))), self._latt.typeIndex)
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

    def test_spglib_input(self):
        _ip = self._latt.get_spglib_input()
        self.assertTupleEqual((self._latt.cell, self._latt.pos, self._latt.typeIndex), _ip)


class lattice_raise(ut.TestCase):
    
    def test_bad_cell(self):
        _cell = [[5.0, 0.0, 0.0],
                 [0.0, 0.0, 5.0]]
        _atoms = ["C"]
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(LatticeError, lattice, _cell, _atoms, _pos)
        _cell = [[5.0, 0.0],
                 [0.0, 5.0, 0.0],
                 [0.0, 0.0, 5.0]]
        self.assertRaises(LatticeError, lattice, _cell, _atoms, _pos)

    def test_bad_atoms_pos(self):
        _cell = [[5.0, 0.0, 0.0],
                 [0.0, 5.0, 0.0],
                 [0.0, 0.0, 5.0]]

        _atoms = []
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(LatticeError, lattice, _cell, _atoms, _pos)

        _atoms = ["C"]
        _pos = [0.0, 0.0, 0.0]
        self.assertRaises(LatticeError, lattice, _cell, _atoms, _pos)

        _atoms = ["C"]
        _pos = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]
        self.assertRaises(LatticeError, lattice, _cell, _atoms, _pos)

        _atoms = ["C", "C"]
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(LatticeError, lattice, _cell, _atoms, _pos)


class lattice_factory_method(ut.TestCase):
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

    def test_read_from_json(self):
        import os
        import tempfile
        import json

        self.assertRaisesRegex(LatticeError, "JSON file not found: None", \
            lattice.read_from_json, None)
        self.assertRaisesRegex(LatticeError, "JSON file not found: /abcdefg.json", \
            lattice.read_from_json, "/abcdefg.json")
        # raise for invalid json
        _tf = tempfile.NamedTemporaryFile()
        with open(_tf.name, 'w') as h:
            json.dump({}, h)
        self.assertRaisesRegex(LatticeError, "invalid JSON file for lattice: {}".format(_tf.name), \
            lattice.read_from_json, _tf.name)

        _dict = {"cell": [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]}
        with open(_tf.name, 'w') as h:
            json.dump(_dict, h)
        self.assertRaisesRegex(LatticeError, "invalid JSON file for lattice: {}. No {}".format(_tf.name, "atoms"), \
            lattice.read_from_json, _tf.name)

        _dict = {
            "cell": [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]], 
            "atoms": ["C"],
            }
        with open(_tf.name, 'w') as h:
            json.dump(_dict, h)
        self.assertRaisesRegex(LatticeError, "invalid JSON file for lattice: {}. No {}".format(_tf.name, "pos"), \
            lattice.read_from_json, _tf.name)
        _tf.close()
        # test one file in testdata, Latt_1.json is tested here
        _path = os.path.join(os.path.dirname(__file__), '..', 'testdata', 'Latt_1.json')
        _latt = lattice.read_from_json(_path)
        self.assertEqual(_latt.unit, "ang")
        self.assertEqual(_latt.coordSys, "D")
        self.assertEqual(_latt.natoms, 2)
        self.assertListEqual(_latt.atoms, ["C", "C"])



class lattice_select_dynamics(ut.TestCase):
    '''Test the functionality of selective dynamics
    '''
    def test_fix_all(self):
        _pc = lattice.bravis_cP("C")
        self.assertFalse(_pc.useSelDyn)
        _pc = lattice.bravis_cP("C", allRelax=False)
        self.assertTrue(_pc.useSelDyn)
        self.assertListEqual([False,]*3, _pc.sdFlags(0))
        _pc = lattice.bravis_cI("C", allRelax=False)
        self.assertTrue(_pc.useSelDyn)
        self.assertListEqual([[False,]*3,]*2, _pc.sdFlags())
    
    def test_fix_some(self):
        _pc = lattice.bravis_cF("C", selectDyn={1:[False, True, True]})
        self.assertListEqual([True,]*3, _pc.sdFlags(0))
        self.assertListEqual([False,True,True,], _pc.sdFlags(1))

    def test_fix_by_set_method(self):
        _pc = lattice.bravis_cF("C")
        _pc.set_fix(0, 1)
        self.assertListEqual([False,False,False], _pc.sdFlags(0))
        self.assertListEqual([False,False,False], _pc.sdFlags(1))
        _pc.set_fix(2, axis=1)
        self.assertListEqual([False,True,True], _pc.sdFlags(2))
        _pc.set_fix(3, axis=[2, 3])
        self.assertListEqual([True,False,False], _pc.sdFlags(3))

    def test_relax_by_set_method(self):
        _pc = lattice.bravis_cF("C", allRelax=False)
        _pc.set_relax(0, 1)
        self.assertListEqual([True,True,True], _pc.sdFlags(0))
        self.assertListEqual([True,True,True], _pc.sdFlags(1))
        _pc.set_relax(2, axis=1)
        self.assertListEqual([True,False,False], _pc.sdFlags(2))
        _pc.set_relax(3, axis=[2, 3])
        self.assertListEqual([False,True,True], _pc.sdFlags(3))

class lattice_sort(ut.TestCase):
    '''Test the sorting functionality of lattice
    '''

    def test_direct_switch_cscl(self):
        _cell = [[1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]]
        _atoms = ["Cl", "Cs"]
        _pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        # Cs atom is fixed
        _fix = [False, False, False]

        _latt = lattice(_cell, _atoms, _pos, selectDyn={1: _fix})
        self.assertListEqual([0], _latt.get_sym_index("Cl"))
        self.assertListEqual([0], _latt["Cl"])
        self.assertListEqual([1], _latt.get_sym_index("Cs"))
        self.assertListEqual([1], _latt["Cs"])
        _latt._switch_two_atom_index(0, 1)
        self.assertListEqual(_latt.atoms, ["Cs", "Cl"])
        self.assertListEqual([0], _latt.get_sym_index("Cs"))
        self.assertListEqual([1], _latt.get_sym_index("Cl"))
        self.assertListEqual(_fix, _latt.sdFlags(0))

    def test_sanitize_atoms_sic(self):
        _cell = [[1.0, 0.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0]]
        _atoms = ["Si", "C", "Si", "Si", "C", "C", "Si", "C"]
        _pos = [[0.0, 0.0, 0.0],     #Si 
                [0.25, 0.25, 0.25],  #C
                [0.0, 0.5, 0.5],     #Si
                [0.5, 0.0, 0.5],     #Si
                [0.25, 0.75, 0.75],  #C
                [0.75, 0.25, 0.75],  #C
                [0.5, 0.5, 0.0],     #Si
                [0.75, 0.75, 0.25]]  #C
        _posSanitied = [[0.0, 0.0, 0.0],     #Si 
                        [0.0, 0.5, 0.5],     #Si
                        [0.5, 0.0, 0.5],     #Si
                        [0.5, 0.5, 0.0],     #Si
                        [0.25, 0.25, 0.25],  #C
                        [0.25, 0.75, 0.75],  #C
                        [0.75, 0.25, 0.75],  #C
                        [0.75, 0.75, 0.25]]  #C
        _latt = lattice(_cell, _atoms, _pos, selectDyn={2:[False, False, False]})
        # _latt._sanitize_atoms()
        self.assertListEqual(list(sorted(_atoms, reverse=True)), _latt.atoms)
        self.assertDictEqual({0: 'Si', 1: 'C'}, _latt.typeMapping)
        self.assertTrue(np.array_equal(_latt.pos, np.array(_posSanitied, dtype=_latt._dtype)))

    def test_sort_pos_sic(self):
        '''Test sorting atoms and their positions in SiC
        '''
        _cell = [[1.0, 0.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0]]
        _atoms = ["Si", "Si", "Si", "Si", "C", "C", "C", "C"]
        _pos = [[0.0, 0.0, 0.0],     #Si 
                [0.0, 0.5, 0.5],     #Si
                [0.5, 0.0, 0.5],     #Si
                [0.5, 0.5, 0.0],     #Si
                [0.25, 0.25, 0.25],  #C
                [0.25, 0.75, 0.75],  #C
                [0.75, 0.25, 0.75],  #C
                [0.75, 0.75, 0.25]]  #C
        _posSorted = [0.5, 0.5, 0.0, 0.0, 0.75, 0.75, 0.25, 0.25]
        _posSortedRev = [0.0, 0.0, 0.5, 0.5, 0.25, 0.25, 0.75, 0.75]
        _latt = lattice(_cell, _atoms, _pos)
        # no need to sanitize atoms
        self.assertListEqual(_atoms, _latt.atoms)
        self.assertDictEqual({0: 'Si', 1: 'C'}, _latt.typeMapping)
        for _axis in range(3):
            _latt.sort_pos(axis=_axis+1)
            self.assertTrue(np.array_equal(np.array(_posSorted, dtype=_latt._dtype), \
                _latt.pos[:,_axis]))
            _latt.sort_pos(axis=_axis+1, reverse=True)
            self.assertTrue(np.array_equal(np.array(_posSortedRev, dtype=_latt._dtype), \
                _latt.pos[:,_axis]))


class test_lattice_utils(ut.TestCase):
    '''Test the utils for ``lattice`` use
    '''
    
    def test_periodic_duplates_in_cell(self):

        _dupcs, _n = periodic_duplicates_in_cell([0.2, 0.4, 0.8])
        self.assertEqual(1, _n)
        self.assertTupleEqual(([0.2, 0.4, 0.8],), _dupcs)
        _dupcs, _n = periodic_duplicates_in_cell([0.2, 0.4, 0])
        self.assertEqual(2, _n)
        self.assertTupleEqual(([0.2, 0.4, 0], [0.2, 0.4, 1.0]), _dupcs)
        _dupcs, _n = periodic_duplicates_in_cell([0, 0.4, 0])
        self.assertEqual(4, _n)
        self.assertTupleEqual(([0, 0.4, 0], [1.0, 0.4, 0], [0, 0.4, 1.0], [1.0, 0.4, 1.0]), _dupcs)
        _dupcs, _n = periodic_duplicates_in_cell([0, 0, 0])
        self.assertEqual(8, _n)
        self.assertTupleEqual(([0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [1.0, 1.0, 0], [0, 0, 1.0], [1.0, 0, 1.0], [0, 1.0, 1.0], [1.0, 1.0, 1.0]), _dupcs)
        self.assertRaises(AssertionError, periodic_duplicates_in_cell, [1.0,0.0,0.0])
        self.assertRaises(AssertionError, periodic_duplicates_in_cell, [1.1,0.0,0.0])
    
    def test_axis_list(self):
        self.assertSetEqual(set([1,2,3]), set(axis_list(0)))
        self.assertSetEqual(set([1,]), set(axis_list(1)))
        self.assertSetEqual(set([1,2]), set(axis_list([1,2])))

    def test_atoms_from_sym_nat(self):
        _atoms = atoms_from_sym_nat(["C", "Al", "F"], [2, 3, 1])
        self.assertListEqual(_atoms, ["C", "C", "Al", "Al", "Al", "F"])

    def test_sym_nat_from_atoms(self):
        _syms, _nats = sym_nat_from_atoms(["C", "Al", "Al", "C", "Al", "F"])
        self.assertListEqual(_syms, ["C", "Al", "F"])
        self.assertListEqual(_nats, [2, 3, 1])


if __name__ == "__main__":
    ut.main()
