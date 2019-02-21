#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import unittest as ut
import logging
from mykit.vasp.poscar import poscar, poscarError
from mykit.core.lattice import lattice

class poscar_build_test(ut.TestCase):
    
    def test_init_from_cell_atoms_pos(self):
        _latt = lattice.bravis_cP("C", aLatt=5.0)
        _poscar = poscar(*_latt.get_latt(), coordSys=_latt.coordSys, unit=_latt.unit)
        _cell = _poscar.get_latt()[0]
        self.assertAlmostEqual(_cell[0,0], 5.0)
        self.assertEqual(len(_poscar), 1)

        _latt = lattice.bravis_cI("C", aLatt=5.0, unit="au")
        _poscar = poscar(*_latt.get_latt(), coordSys=_latt.coordSys, unit=_latt.unit)
        self.assertEqual(_poscar.unit, "au")
        self.assertEqual(len(_poscar), 2)

        _latt = lattice.bravis_cF("C", aLatt=5.0, coordSys="C")
        _poscar = poscar(*_latt.get_latt(), coordSys=_latt.coordSys, unit=_latt.unit)
        self.assertEqual(_poscar.coordSys, "C")
        self.assertEqual(len(_poscar), 4)

    def test_init_from_lattice(self):
        _latt = lattice.bravis_cP("C", aLatt=5.0)
        _poscar = poscar.create_from_lattice(_latt)
        self.assertAlmostEqual(_poscar.get_latt()[0][0,0], 5.0)
        self.assertEqual(len(_poscar), 1)

        _latt = lattice.bravis_cI("C", aLatt=5.0, unit="au")
        _poscar = poscar.create_from_lattice(_latt)
        self.assertEqual(_poscar.unit, "au")
        self.assertEqual(len(_poscar), 2)
        
        _latt = lattice.bravis_cF("C", aLatt=5.0, coordSys="C")
        _poscar = poscar.create_from_lattice(_latt)
        self.assertEqual(_poscar.coordSys, "C")
        self.assertEqual(len(_poscar), 4)

    def test_init_from_file(self):

        _countGood = 0
        _countBad = 0
        __poscarDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../testdata')
        if os.path.isdir(__poscarDir):
            for _f in os.listdir(__poscarDir):
                if re.match('^POSCAR_*', _f):
                    # TODO verify the reads
                    _countGood += 1
                    _i = _f.split('_')[1]
                    _path = os.path.join(__poscarDir, _f)
                    _poscar = poscar.read_from_file(_path)
                if re.match('^bad_POSCAR_*', _f):
                    _countBad += 1
                    _i = _f.split('_')[2]
                    _path = os.path.join(__poscarDir, _f)
                    self.assertRaises(poscarError, poscar.read_from_file, _path)
        print("{} good POSCARs readed, {} bad POSCARs readed".format(_countGood, _countBad))


if __name__ == "__main__":
    ut.main()