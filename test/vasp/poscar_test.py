#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import logging
import os
import re
import unittest as ut

import numpy as np

from mykit.core.lattice import lattice
from mykit.vasp.poscar import PoscarError, poscar


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

    def test_read_from_file(self):

        _countGood = 0
        _countBad = 0
        _countVerified = 0
        __poscarDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../testdata')
        if os.path.isdir(__poscarDir):
            for _f in os.listdir(__poscarDir):
                if re.match('^POSCAR_[0-9]+$', _f):
                    _countGood += 1
                    _i = _f.split('_')[1]
                    _path = os.path.join(__poscarDir, _f)
                    _poscar = poscar.read_from_file(_path)
                    _verifyJson = os.path.join(__poscarDir, 'Latt_'+_i+'.json')
                    if os.path.isfile(_verifyJson):
                        _vs = _verify_poscar_by_json(self, _poscar, _verifyJson)
                        if _vs:
                            _countVerified += 1
                if re.match('^bad_POSCAR_[0-9]+$', _f):
                    _countBad += 1
                    _i = _f.split('_')[2]
                    _path = os.path.join(__poscarDir, _f)
                    self.assertRaises(PoscarError, poscar.read_from_file, _path)
        print("{} good POSCARs readed ({} verified by JSON file). {} bad POSCARs raised.".format(_countGood, _countVerified, _countBad))


def _verify_poscar_by_json(tc, pc, pathJson):
    '''
    '''
    assert isinstance(tc, ut.TestCase)
    assert isinstance(pc, poscar)
    with open(pathJson, 'r') as f:
        _vDict = json.load(f)

    _cellPos, _atomsPos, _posPos = pc.get_latt()
    _kwargsPos = pc.get_kwargs()
    verifyMsg = "POSCAR verification failed: {}".format(pathJson)
    if "cell" in _vDict:
        tc.assertTrue(np.array_equal(pc.cell, np.array(_vDict["cell"], dtype=pc._dtype)), msg=verifyMsg)
    if "sdFlags" in _vDict:
        _sdfDict = _vDict["sdFlags"]
        # print(_sdfDict.keys())
        for _a in _sdfDict:
            _ind = int(_a)
            tc.assertListEqual(_sdfDict[_a], pc.sdFlags(_ind), msg=verifyMsg+", sdFlag: atom {}".format(_a))
    if "coordSys" in _vDict:
        tc.assertEqual(_vDict["coordSys"], pc.coordSys, msg=verifyMsg)
    # atoms after arrangement sanitization
    if "atoms" in _vDict:
        tc.assertListEqual(_vDict["atoms"], pc.atoms, msg=verifyMsg)
    if "comment" in _vDict:
        tc.assertEqual(_vDict["comment"], pc.comment)
    if "pos" in _vDict:
        print(pc.pos)
        print(_vDict["pos"])
        tc.assertTrue(np.array_equal(pc.pos, np.array(_vDict["pos"], dtype=pc._dtype)))
    # put sort check at the end
    if "sorted" in _vDict:
        for _axis in range(3):
            pc.sort_pos(axis=_axis+1)
            tc.assertTrue(np.array_equal(pc.pos[:, _axis], \
                np.array(_vDict["sorted"][_axis], dtype=pc._dtype)))
    if "sortedRev" in _vDict:
        for _axis in range(3):
            pc.sort_pos(axis=_axis+1, reverse=True)
            tc.assertTrue(np.array_equal(pc.pos[:, _axis], \
                np.array(_vDict["sortedRev"][_axis], dtype=pc._dtype)))
    return True


if __name__ == "__main__":
    ut.main()
