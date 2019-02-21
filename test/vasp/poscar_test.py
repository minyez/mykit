#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut
from mykit.vasp.poscar import poscar
from mykit.core.lattice import lattice

class build_test(ut.TestCase):
    
    def test_init_from_cell_atoms_pos(self):
        _latt = lattice.bravis_cP("C", aLatt=5.0)
        _poscar = poscar(_latt.cell, _latt.atoms, _latt.pos, coordSys=_latt.coordSys, unit=_latt.unit)
        self.assertAlmostEqual(_poscar.cell[0,0], 5.0)
        self.assertEqual(len(_poscar), 1)
        _latt = lattice.bravis_cI("C", aLatt=5.0, unit="au")
        _poscar = poscar(_latt.cell, _latt.atoms, _latt.pos, coordSys=_latt.coordSys, unit=_latt.unit)
        self.assertEqual(_poscar.unit, "au")
        self.assertEqual(len(_poscar), 2)
        _latt = lattice.bravis_cF("C", aLatt=5.0, coordSys="C")
        _poscar = poscar(_latt.cell, _latt.atoms, _latt.pos, coordSys=_latt.coordSys, unit=_latt.unit)
        self.assertEqual(_poscar.coordSys, "C")
        self.assertEqual(len(_poscar), 4)

    def test_init_from_lattice(self):
        _latt = lattice.bravis_cP("C", aLatt=5.0)
        _poscar = poscar.create_from_lattice(_latt)
        self.assertAlmostEqual(_poscar.cell[0,0], 5.0)
        self.assertEqual(len(_poscar), 1)
        


if __name__ == "__main__":
    ut.main()