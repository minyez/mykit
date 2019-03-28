#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import unittest as ut

import numpy as np

from mykit.core.cell import (Cell, CellError, atoms_from_sym_nat, axis_list,
                             periodic_duplicates_in_cell, sym_nat_from_atoms)
from mykit.core.constants import ANG2AU, AU2ANG, PI
from mykit.core.utils import get_dirpath, get_matched_files


class simple_cubic_lattice(ut.TestCase):

    _a = 5.0
    _latt = [[_a, 0.0, 0.0],
             [0.0, _a, 0.0],
             [0.0, 0.0, _a]]
    _atoms = ["C"]
    _frac = 0.5
    _pos = [[_frac, 0.0, 0.0]]

    _cell = Cell(_latt, _atoms, _pos, unit="ang", coordSys="D")
    
    def test_properties(self):
        # unit
        self.assertEqual(self._cell.unit, "ang")
        # coordinate system
        self.assertEqual(self._cell.coordSys, "D")
        # real space vectors
        self.assertTrue(np.array_equal(self._cell.a, self._latt))
        self.assertAlmostEqual(pow(self._a, 3), self._cell.vol)
        self.assertTrue(np.array_equal(self._cell.center, [0.5,0.5,0.5]))
        lattConsts = self._cell.lattConsts
        self.assertTupleEqual(lattConsts, \
            (self._a, self._a, self._a, 90.0, 90.0, 90.0))
        # Reciprocal
        recpIn2Pi = np.array(self._latt, dtype=self._cell._dtype)/self._a**2
        self.assertTrue(np.allclose(self._cell.bIn2Pi, recpIn2Pi))
        self.assertTrue(np.allclose(self._cell.b, recpIn2Pi * 2.0 * PI))
        self.assertTrue(np.allclose(self._cell.blen, (2.0*PI/self._a,)*3))
        # atom types
        self.assertEqual(1, self._cell.natoms)
        self.assertListEqual(["C",], self._cell.atomTypes)
        self.assertDictEqual({0: "C"}, self._cell.typeMapping)
        self.assertListEqual([0,], self._cell.typeIndex)
        
    def test_magic(self):
        self.assertEqual(1, len(self._cell))
        self.assertTupleEqual(tuple(self._pos[0]), tuple(self._cell[0]))

    def test_unit_conv(self):
        # ang2au
        self._cell.unit = 'au'
        _latt = self._cell.get_cell()[0]
        if self._cell._dtype == 'float32':
            self.assertAlmostEqual(_latt[0,0], self._a * ANG2AU, places=5)
        else:
            self.assertAlmostEqual(_latt[0,0], self._a * ANG2AU)
        # au2ang
        self._cell.unit = 'ang'
        _latt = self._cell.get_cell()[0]
        self.assertEqual(_latt[0,0], self._a)

    def test_coord_conv(self):
        # direct2cart
        self._cell.coordSys = 'C'
        self.assertTupleEqual(tuple(self._cell[0]), (self._frac * self._a, 0.0, 0.0))
        self.assertEqual("C", self._cell.coordSys)
        # cart2direct
        self._cell.coordSys = 'D'
        self.assertEqual(self._cell[0][0], self._frac)
        self.assertEqual("D", self._cell.coordSys)

    def test_spglib_input(self):
        _ip = self._cell.get_spglib_input()
        self.assertTupleEqual((self._cell.latt, self._cell.pos, self._cell.typeIndex), _ip)


class cell_raise(ut.TestCase):
    
    def test_bad_cell(self):
        _latt = [[5.0, 0.0, 0.0],
                 [0.0, 0.0, 5.0]]
        _atoms = ["C"]
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(CellError, Cell, _latt, _atoms, _pos)
        _latt = [[5.0, 0.0],
                 [0.0, 5.0, 0.0],
                 [0.0, 0.0, 5.0]]
        self.assertRaises(CellError, Cell, _latt, _atoms, _pos)

    def test_bad_atoms_pos(self):
        _latt = [[5.0, 0.0, 0.0],
                 [0.0, 5.0, 0.0],
                 [0.0, 0.0, 5.0]]

        _atoms = []
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(CellError, Cell, _latt, _atoms, _pos)

        _atoms = ["C"]
        _pos = [0.0, 0.0, 0.0]
        self.assertRaises(CellError, Cell, _latt, _atoms, _pos)

        _atoms = ["C"]
        _pos = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]
        self.assertRaises(CellError, Cell, _latt, _atoms, _pos)

        _atoms = ["C", "C"]
        _pos = [[0.0, 0.0, 0.0]]
        self.assertRaises(CellError, Cell, _latt, _atoms, _pos)


class cell_factory_method(ut.TestCase):
    '''Test the class methods to generate commonly used lattice structure'''

    def test_bravais_cubic(self):
        _pc = Cell.bravais_cP("C", a=5.0, coordSys="D")
        self.assertEqual(1, len(_pc))
        self.assertEqual("D", _pc.coordSys)
        _bcc = Cell.bravais_cI("C", a=5.0, primitive=False, unit="au")
        self.assertEqual("au", _bcc.unit)
        self.assertEqual(2, len(_bcc))
        _fcc = Cell.bravais_cF("C", a=5.0, primitive=False)
        self.assertEqual(4, len(_fcc))
        # primitive cell
        _pbcc = Cell.bravais_cI("C", a=5.0, primitive=True)
        self.assertEqual(1, len(_pbcc))
        self.assertAlmostEqual(5.0*np.sqrt(3.0)/2.0, _pbcc.alen[0])
        _pfcc = Cell.bravais_cF("C", a=5.0, primitive=True)
        self.assertEqual(1, len(_pfcc))
        self.assertAlmostEqual(5.0*np.sqrt(0.5), _pfcc.alen[0])

    def test_bravais_orth(self):
        oP = Cell.bravais_oP("C", a=1.0, b=2.0, c=3.0)
        self.assertEqual(len(oP), 1)
        self.assertEqual(oP.vol, 6.0)
        oI = Cell.bravais_oI("C")
        self.assertEqual(len(oI), 2)
        oF = Cell.bravais_oF("C")
        self.assertEqual(len(oF), 4)

    def test_typical_sysmtes(self):
        # both conventional and primitive
        for p in [True, False]:
            Cell.diamond("C", primitive=p)
            Cell.anatase("Ti", "O", primitive=p)
            Cell.rutile("Ti", "O", primitive=p)
            Cell.zincblende("Zn", "O", primitive=p)
        # primitive only
        Cell.perovskite("Ca", "Ti", "O")
        Cell.wurtzite("Zn", "O")
        Cell.pyrite()
        Cell.marcasite()

    def test_read_from_json(self):
        import os
        import tempfile
        import json

        self.assertRaisesRegex(CellError, "JSON file not found: None", \
            Cell.read_from_json, None)
        self.assertRaisesRegex(CellError, "JSON file not found: /abcdefg.json", \
            Cell.read_from_json, "/abcdefg.json")
        # raise for invalid json
        _tf = tempfile.NamedTemporaryFile()
        with open(_tf.name, 'w') as h:
            json.dump({}, h)
        self.assertRaisesRegex(CellError, "invalid JSON file for cell: {}".format(_tf.name), \
            Cell.read_from_json, _tf.name)

        _dict = {"latt": [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]}
        with open(_tf.name, 'w') as h:
            json.dump(_dict, h)
        self.assertRaisesRegex(CellError, "invalid JSON file for cell: {}. No {}".format(_tf.name, "atoms"), \
            Cell.read_from_json, _tf.name)

        _dict = {
            "latt": [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]], 
            "atoms": ["C"],
            }
        with open(_tf.name, 'w') as h:
            json.dump(_dict, h)
        self.assertRaisesRegex(CellError, "invalid JSON file for cell: {}. No {}".format(_tf.name, "pos"), \
            Cell.read_from_json, _tf.name)

        # JSON with factory key
        _dict = {
            "factory": "zincblende", 
            "atom1": "Zn", 
            "a": 8.0, "unit": "au",
            }
        with open(_tf.name, 'w') as h:
            json.dump(_dict, h)
        self.assertRaisesRegex(CellError, "Required key not found in JSON: atom2", \
            Cell.read_from_json, _tf.name)
        # add atom2 and dump again
        _dict["atom2"] = "O"
        with open(_tf.name, 'w') as h:
            json.dump(_dict, h)
        _cell = Cell.read_from_json(_tf.name)
        self.assertEqual(_cell.unit, 'au')
        self.assertEqual(_cell.comment, "Zincblende ZnO")
        self.assertAlmostEqual(512, _cell.vol)
        _tf.close()

        # test one file in testdata, Cell_1.json is tested here
        _path = os.path.join(os.path.dirname(__file__), '..', 'testdata', 'Cell_1.json')
        _cell = Cell.read_from_json(_path)
        self.assertEqual(_cell.unit, "ang")
        self.assertEqual(_cell.coordSys, "D")
        self.assertEqual(_cell.natoms, 2)
        self.assertListEqual(_cell.atoms, ["C", "C"])
    
    def test_read_from_cif(self):
        dataDir = os.path.join(get_dirpath(__file__), '..', 'testdata')
        for fn in os.listdir(dataDir):
            if fn.endswith('.cif'):
                cif = os.path.join(dataDir, fn)
                Cell.read_from_cif(cif)


class cell_select_dynamics(ut.TestCase):
    '''Test the functionality of selective dynamics
    '''
    def test_fix_all(self):
        _c = Cell.bravais_cP("C")
        self.assertFalse(_c.useSelDyn)
        _c = Cell.bravais_cP("C", allRelax=False)
        self.assertTrue(_c.useSelDyn)
        self.assertListEqual([False,]*3, _c.sdFlags(0))
        _c = Cell.bravais_cI("C", allRelax=False, primitive=False)
        self.assertTrue(_c.useSelDyn)
        self.assertListEqual([[False,]*3,]*2, _c.sdFlags())

    def test_relax_all(self):
        _c = Cell.bravais_cI("C", allRelax=False, primitive=False, \
            selectDyn={1: [True, False, True]})
        _c.relax_all()
        self.assertListEqual(_c.sdFlags(1), [True, True, True])
        self.assertFalse(_c.useSelDyn)
    
    def test_fix_some(self):
        _pc = Cell.bravais_cF("C", selectDyn={1:[False, True, True]})
        self.assertListEqual([True,]*3, _pc.sdFlags(0))
        self.assertListEqual([False,True,True,], _pc.sdFlags(1))

    def test_fix_by_set_method(self):
        _pc = Cell.bravais_cF("C")
        _pc.set_fix(0, 1)
        self.assertListEqual([False,False,False], _pc.sdFlags(0))
        self.assertListEqual([False,False,False], _pc.sdFlags(1))
        _pc.set_fix(2, axis=1)
        self.assertListEqual([False,True,True], _pc.sdFlags(2))
        _pc.set_fix(3, axis=[2, 3])
        self.assertListEqual([True,False,False], _pc.sdFlags(3))

    def test_relax_by_set_method(self):
        _pc = Cell.bravais_cF("C", allRelax=False)
        _pc.set_relax(0, 1)
        self.assertListEqual([True,True,True], _pc.sdFlags(0))
        self.assertListEqual([True,True,True], _pc.sdFlags(1))
        _pc.set_relax(2, axis=1)
        self.assertListEqual([True,False,False], _pc.sdFlags(2))
        _pc.set_relax(3, axis=[2, 3])
        self.assertListEqual([False,True,True], _pc.sdFlags(3))

class cell_sort(ut.TestCase):
    '''Test the sorting functionality of Cell
    '''

    def test_direct_switch_cscl(self):
        _latt = [[1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]]
        _atoms = ["Cl", "Cs"]
        _pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        # Cs atom is fixed
        _fix = [False, False, False]

        _cell = Cell(_latt, _atoms, _pos, selectDyn={1: _fix})
        self.assertListEqual([0], _cell.get_sym_index("Cl"))
        self.assertListEqual([0], _cell["Cl"])
        self.assertListEqual([1], _cell.get_sym_index("Cs"))
        self.assertListEqual([1], _cell["Cs"])
        _cell._switch_two_atom_index(0, 1)
        self.assertListEqual(_cell.atoms, ["Cs", "Cl"])
        self.assertListEqual([0], _cell.get_sym_index("Cs"))
        self.assertListEqual([1], _cell.get_sym_index("Cl"))
        self.assertListEqual(_fix, _cell.sdFlags(0))

    def test_sanitize_atoms_sic(self):
        _latt = [[1.0, 0.0, 0.0],
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
        SiC = Cell(_latt, _atoms, _pos, \
            selectDyn={2:[False, False, False]})
        # _latt._sanitize_atoms()
        self.assertListEqual(list(sorted(_atoms, reverse=True)), \
            SiC.atoms)
        self.assertDictEqual({0: 'Si', 1: 'C'}, SiC.typeMapping)
        self.assertTrue(np.array_equal(SiC.pos, \
            np.array(_posSanitied, dtype=SiC._dtype)))
        self.assertListEqual([False, False, False], SiC.sdFlags(1))

    def test_sort_pos_sic(self):
        '''Test sorting atoms and their positions in SiC
        '''
        _latt = [[1.0, 0.0, 0.0],
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
        SiC = Cell(_latt, _atoms, _pos)
        # no need to sanitize atoms
        self.assertListEqual(_atoms, SiC.atoms)
        self.assertDictEqual({0: 'Si', 1: 'C'}, SiC.typeMapping)
        for _axis in range(3):
            SiC.sort_pos(axis=_axis+1)
            self.assertTrue(np.array_equal(np.array(_posSorted, dtype=SiC._dtype), \
                SiC.pos[:,_axis]))
            SiC.sort_pos(axis=_axis+1, reverse=True)
            self.assertTrue(np.array_equal(np.array(_posSortedRev, dtype=SiC._dtype), \
                SiC.pos[:,_axis]))


class test_cell_manipulation(ut.TestCase):
    '''Test manipulation methods for lattice and atoms
    '''

    def test_add_atom_on_graphene(self):
        '''Test adding atoms in a graphene cell
        '''
        a = 5.2
        _latt = [[a/2, 0.0, 0.0],
                 [-a/2, a/2*np.sqrt(3.0), 0.0],
                 [0.0, 0.0, 15.0]]
        _atoms = ["C", "C",]
        _pos = [[0.0, 0.0, 0.5],     
                [1.0/3,2.0/3,0.5],]
        gp = Cell(_latt, _atoms, _pos)
        self.assertRaisesRegex(CellError, \
            r"Invalid coordinate: *", \
            gp.add_atom, "H", [0.2, 0.3])
        self.assertRaisesRegex(CellError, \
            "atom should be string, received <class 'int'>", \
            gp.add_atom, 1, [0.2,0.3,0.4])
        gp.fix_all()
        gp.add_atom("H", [0.0,0.0,0.6], sdFlag=[False, False, True])
        self.assertEqual(gp.natoms, 3)
        self.assertListEqual(gp.atoms, ['C','C','H'])
        self.assertDictEqual(gp.typeMapping, {0: 'C', 1: 'H'})
        self.assertListEqual(gp.sdFlags(2), [False, False, True])

    def test_atom_arrange_after_add_atom(self):
        '''Test if the atoms are correctly rearranged
        after adding new atom
        '''
        a = 2.0
        _latt = [[a, 0.0, 0.0],
                 [0.0, a, 0.0],
                 [0.0, 0.0, a]]
        _atoms = ["Na", "Cl",]
        _pos = [[0.0, 0.0, 0.0],     
                [0.5, 0.5, 0.5],]
        brokenNaCl = Cell(_latt, _atoms, _pos)
        brokenNaCl.add_atom('Na', [0.0,0.5,0.5])
        self.assertListEqual(brokenNaCl.atoms, ['Na','Na','Cl'])
        brokenNaCl.add_atom('Na', [0.5,0.0,0.5])
        self.assertListEqual(brokenNaCl.atoms, ['Na','Na','Na','Cl'])
        brokenNaCl.add_atom('Cl', [0.5,0.0,0.0])
        self.assertListEqual(brokenNaCl.atoms, ['Na','Na','Na','Cl','Cl'])
        self.assertTrue(np.array_equal(brokenNaCl.pos, \
            np.array([[0.0,0.0,0.0],
                      [0.0,0.5,0.5],
                      [0.5,0.0,0.5],
                      [0.5,0.5,0.5],
                      [0.5,0.0,0.0],], dtype=brokenNaCl._dtype)))


class test_cell_utils(ut.TestCase):
    '''Test the utility functions for ``Cell`` use
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
