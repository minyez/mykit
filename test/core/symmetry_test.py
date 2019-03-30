#!/usr/bin/env python3
# coding = utf-8

import re
import unittest as ut

import numpy as np

from mykit.core.cell import Cell
from mykit.core.metadata._spk import (_special_kpoints, cond_a_lt_b,
                                      cond_abc_invsq, cond_any, cond_c_lt_a,
                                      cond_curt_a_lt_sqrt_c,
                                      cond_max_bc_noneq_left)
from mykit.core.symmetry import (SpaceGroup, SpecialKpoints, SymmetryError,
                                 _check_valid_custom_ksym_dict,
                                 _check_valid_spg_id,
                                 _spglib_check_cell_and_coordSys, get_ibzkpt,
                                 get_inequiv, get_spacegroup,
                                 get_sym_operations, is_primitive, primitize,
                                 standardize)

NKSETS_COND = {
    cond_any: 1,
    cond_a_lt_b: 2,
    cond_c_lt_a: 2,
    cond_curt_a_lt_sqrt_c: 2,
    cond_max_bc_noneq_left: 2,
    cond_abc_invsq: 3,
}


class test_crystal_symmetries(ut.TestCase):
    '''
    '''
    # rutile TiO2
    rTiO2 = Cell.rutile("Ti", "O")
    # w-ZnS
    wZnS = Cell.wurtzite("Zn", "S")
    # zincblende ZnO
    zbPrim = Cell.zincblende("Zn", "O", primitive=True)
    zbConv = Cell.zincblende("Zn", "O")
    # diamond
    dcPrim = Cell.diamond("C", aLatt=3.25, primitive=True)
    dcConv = Cell.diamond("C", aLatt=3.25, primitive=False)

    def test_sym_operations(self):
        self.assertEqual(16, len(get_sym_operations(self.rTiO2)["rotations"]))
        self.assertEqual(16, len(get_sym_operations(self.rTiO2)["translations"]))

    def test_is_primitive(self):
        self.assertTrue(is_primitive(self.zbPrim))
        self.assertTrue(is_primitive(self.rTiO2))
        self.assertTrue(is_primitive(self.wZnS))
        self.assertTrue(is_primitive(self.dcPrim))
        self.assertFalse(is_primitive(self.zbConv))
        self.assertFalse(is_primitive(self.dcConv))

    def test_well_known_spg_of_sysmtems(self):
        self.assertTupleEqual((136, 'P4_2/mnm'), get_spacegroup(self.rTiO2))
        self.assertTupleEqual((186, 'P6_3mc'), get_spacegroup(self.wZnS))
        self.assertTupleEqual((216, 'F-43m'), get_spacegroup(self.zbPrim))
        self.assertTupleEqual((227, 'Fd-3m'), get_spacegroup(self.dcPrim))
    
    def test_get_inequiv(self):
        a, _i = get_inequiv(self.dcPrim)
        self.assertListEqual(['C', ], a)
        a, _i = get_inequiv(self.zbPrim)
        self.assertListEqual(['Zn', 'O'], a)

    def test_zincblende_prim(self):
        # irreducible kpoints
        get_ibzkpt(self.zbConv, [2, 2, 2])
        # primitize
        zbPrim = primitize(self.zbConv)
        self.assertIsInstance(zbPrim, Cell)
        self.assertTrue(is_primitive(zbPrim))
        # standardize
        std = standardize(self.zbConv, primitive=False)
        self.assertIsInstance(std, Cell)
        # for zincblende, the standard cell is conventional
        self.assertFalse(is_primitive(std))
        std = standardize(self.zbConv, primitive=True)
        self.assertIsInstance(std, Cell)
        self.assertTrue(is_primitive(std))


class test_space_group(ut.TestCase):

    testPair = [
        (1, 'P1'), (22, 'F222'), (64, 'Cmce'), (123, 'P4mmm'),
        (161, 'R3c'), (168, 'P6_3mc'), (208, 'P4232'), (230, 'Ia-3d')
    ]
    
    def test_symbols(self):
        self.assertEqual(len(SpaceGroup.symbols), 230)

    def test_get_spg_symbol(self):
        for i, symbol in self.testPair:
            self.assertEqual(SpaceGroup.get_spg_symbol(i), symbol)
        # test raises
        self.assertRaisesRegex(SymmetryError, r"string received. id should be int",
                               SpaceGroup.get_spg_symbol, 'P1')
        self.assertRaisesRegex(SymmetryError, r"Invalid space group id \(1~230\): 231",
                               SpaceGroup.get_spg_symbol, 231)
        self.assertRaisesRegex(SymmetryError, r"Invalid space group id \(1~230\): -1",
                               SpaceGroup.get_spg_symbol, -1)

    def test_get_spg_id(self):
        for i, symbol in self.testPair:
            self.assertEqual(SpaceGroup.get_spg_id(symbol), i)
        self.assertRaisesRegex(SymmetryError,
                               r"int received. symbol should be string",
                               SpaceGroup.get_spg_id, 1)
        self.assertRaisesRegex(SymmetryError,
                               r"Space group symbol is not found in ITA: Pxxx",
                               SpaceGroup.get_spg_id, 'Pxxx')


class test_spk_dict(ut.TestCase):
    '''Test the integrity of _special_kpoints dictionary
    '''

    def test_cond_func(self):
        info = r"Fail to determine special kpoints set. May need to standardize first."
        self.assertEqual(cond_any(1, 2, 3), 0)

        self.assertEqual(cond_a_lt_b(2, 1, 3), 0)
        self.assertEqual(cond_a_lt_b(1, 2, 3), 1)
        self.assertRaisesRegex(ValueError, info, cond_a_lt_b, 1, 1, 3)

        self.assertEqual(cond_c_lt_a(2, 1, 3), 0)
        self.assertEqual(cond_c_lt_a(3, 1, 2), 1)
        self.assertRaisesRegex(ValueError, info, cond_c_lt_a, 1, 3, 1)

        self.assertEqual(cond_curt_a_lt_sqrt_c(8, 1, 3), 0)
        self.assertEqual(cond_curt_a_lt_sqrt_c(8, 1, 6), 1)
        self.assertRaisesRegex(
            ValueError, info, cond_curt_a_lt_sqrt_c, 8, 1, 4)

        self.assertEqual(cond_max_bc_noneq_left(2, 5, 3), 0)
        self.assertEqual(cond_max_bc_noneq_left(3, 5, 2), 0)
        self.assertEqual(cond_max_bc_noneq_left(2, 3, 5), 1)
        self.assertEqual(cond_max_bc_noneq_left(3, 2, 5), 1)
        self.assertRaisesRegex(
            ValueError, info, cond_max_bc_noneq_left, 2, 5, 2)
        self.assertRaisesRegex(
            ValueError, info, cond_max_bc_noneq_left, 2, 2, 5)
        self.assertRaisesRegex(
            ValueError, info, cond_max_bc_noneq_left, 5, 2, 2)

        self.assertEqual(cond_abc_invsq(0.1, 0.5, 0.5), 0)
        self.assertEqual(cond_abc_invsq(0.5, 0.5, 0.1), 1)
        self.assertEqual(cond_abc_invsq(0.2, 0.2, 0.5), 2)
        self.assertRaisesRegex(ValueError, info, cond_abc_invsq, 0.5, 0.1, 0.5)

    def test_number_of_spg(self):
        for i in range(1, 231):
            self.assertIn(i, _special_kpoints.keys())
        self.assertEqual(230, len(_special_kpoints.keys()))

    def test_has_cond_and_spPrim(self):
        for k, v in _special_kpoints.items():
            self.assertIn("cond", v,
                          msg="missing cond key for spg {}".format(k))
            self.assertIn("spPrim", v,
                          msg="missing spPrim key for spg {}".format(k))
            self.assertTrue(isinstance(v["spPrim"], list),
                            msg="bad spPrim value for spg {}".format(k))

    def test_number_of_sets(self):
        for k, v in _special_kpoints.items():
            self.assertEqual(NKSETS_COND[v["cond"]], len(v["spPrim"]),
                             msg="inconsistent cond and spPrim for spg {}".format(k))

    def test_kpath_value(self):
        '''Check the value of kpath key

        Explicitly, if the ``kpath`` key exists, 

        - its value should be a list of lists with the same len as value of ``spPrim`` 
            key

        - The lists within are lists of strings, each of which is a kpath and should match
            the ``KPATH_PATTERN`` defined in ``kmesh`` module

        - the symbols obtained from decoding the kpath string, should belong to the
            corresponding special kpoints set.
        '''
        from mykit.core.kmesh import KPATH_PATTERN, kpath_decoder
        for k, v in _special_kpoints.items():
            msg = "bad kpath value for spg {}".format(k)
            kpathSets = v.get('kpath', None)
            if kpathSets != None:
                self.assertTrue(isinstance(kpathSets, list))
                nksets = NKSETS_COND[v["cond"]]
                self.assertEqual(nksets, len(kpathSets), msg=msg)
                for i, kpaths in enumerate(kpathSets):
                    if kpaths != []:
                        self.assertTrue(isinstance(kpaths, list), msg=msg)
                        ksyms = v["spPrim"][i].keys()
                        for kpath in kpaths:
                            self.assertTrue(isinstance(kpath, str), msg=msg)
                            for seg in kpath.split():
                                self.assertRegex(seg, KPATH_PATTERN, msg=msg)
                            decoded = kpath_decoder(kpath)
                            for ksym in decoded:
                                self.assertIn(ksym, ksyms, msg=msg)


class test_special_kpoints(ut.TestCase):

    def test_direct_initialize(self):
        '''Test initialize directly or from Cell or Symmetry instance
        '''
        # check P1
        spk = SpecialKpoints(1, (1.0, 1.1, 1.2), True)
        self.assertEqual(1, len(spk.spkCoord))
        self.assertIn("GM", spk.spkCoord)

    def test_convert_kpath(self):
        # check space group 227, primitive
        spkpts = SpecialKpoints(227, (5.0, 5.0, 5.0), True)
        pathGM2X2L = spkpts.convert_kpath("GM-X-L")
        self.assertIn("symbols", pathGM2X2L)
        self.assertIn("coordinates", pathGM2X2L)
        syms = pathGM2X2L["symbols"]
        coords = pathGM2X2L["coordinates"]
        self.assertListEqual(["GM", "X", "X", "L"], syms)
        self.assertEqual(len(coords), 4)
        self.assertTrue(np.array_equal(coords,
                                       np.array([[0.0, 0.0, 0.0],
                                                 [0.5, 0.0, 0.5],
                                                 [0.5, 0.0, 0.5],
                                                 [0.5, 0.5, 0.5]])))
        spkpts = SpecialKpoints(5, (5.0, 5.0, 5.0), False)
        pathGM2Y = spkpts.convert_kpath("GM-Y")
        # Y = (0.5,0.5,0) in prim, (0,1,0)
        coords = pathGM2Y["coordinates"]
        self.assertTrue(np.array_equal(coords,
                                       np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])))

    def test_kpath_facotry_methods(self):
        # Zincblende (spg 216)
        zb = Cell.zincblende("Zn", "O", a=5.0, primitive=True)
        kp = SpecialKpoints.get_kpaths_predef_from_cell(zb)
        self.assertTrue(isinstance(kp, list))
        kp = SpecialKpoints.get_kpaths_predef_from_cell(zb, 0)
        self.assertTrue(isinstance(kp, dict))
        kp = SpecialKpoints.get_kpath_from_cell("GM-L-X", zb)
        symbols = kp["symbols"]
        coords = kp["coordinates"]
        self.assertListEqual(["GM", "L", "L", "X"], symbols)
        self.assertTrue(np.array_equal(np.array(coords),
                                       np.array([[0.0, 0.0, 0.0],
                                                 [0.5, 0.5, 0.5],
                                                 [0.5, 0.5, 0.5],
                                                 [0.5, 0.0, 0.5]])))

    def test_magic(self):
        spkpts = SpecialKpoints(199, (5.0, 5.0, 5.0), True)
        self.assertListEqual([0.0, 0.0, 0.0], spkpts["GM"])
        self.assertListEqual([0.0, 0.0, 0.5], spkpts["N"])
        spkpts["P"] = [0.2, 0.2, 0.2]
        self.assertDictEqual(spkpts._custom, {"P": [0.2, 0.2, 0.2]})

    def test_init_from_cell(self):
        aTiO2Prim = Cell.anatase("Ti", "O", primitive=True)
        # print(aTiO2Prim)
        spk = SpecialKpoints.from_cell(aTiO2Prim)
        self.assertEqual(spk.spgId, 141)
        aTiO2Conv = Cell.anatase("Ti", "O", primitive=False)
        spk = SpecialKpoints.from_cell(aTiO2Conv)
        self.assertEqual(spk.spgId, 141)


class test_symmetry_utils(ut.TestCase):

    def test_check_cell_and_coordSys(self):
        self.assertRaisesRegex(SymmetryError,
                               r"The input should be instance of Cell *",
                               _spglib_check_cell_and_coordSys, 1)
        cF = Cell.bravais_cF("C", a=5.0)
        cF.coordSys = "C"
        # print(cF)
        self.assertRaisesRegex(SymmetryError,
                               r"The coordinate system should be direct. Cartisian found.",
                               _spglib_check_cell_and_coordSys, cF)

    def test_check_valid_spg_id(self):
        self.assertRaisesRegex(SymmetryError,
                               r"string received. id should be int.",
                               _check_valid_spg_id, "1")
        for badId in [1.0, (1,)]:
            self.assertRaisesRegex(SymmetryError,
                                   r"id should be int.",
                                   _check_valid_spg_id, badId)
        for badId in [-10, 235, 400]:
            self.assertRaisesRegex(SymmetryError,
                                   r"Invalid space group id *",
                                   _check_valid_spg_id, badId)

    def test_check_valid_custom_ksym_dict(self):
        self.assertDictEqual({}, _check_valid_custom_ksym_dict(None))
        self.assertRaisesRegex(TypeError,
                               r"custom_symbols must be dict, *",
                               _check_valid_custom_ksym_dict, 1)
        ksymDict = {"P": [0.2, 0.2, 0.2]}
        self.assertDictEqual(ksymDict, _check_valid_custom_ksym_dict(ksymDict))


if __name__ == '__main__':
    ut.main()
