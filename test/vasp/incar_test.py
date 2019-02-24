#!/usr/bin/env python3
# coding = utf-8

import unittest as ut
from mykit.vasp.incar import incar, incarError

class test_direct_set(ut.TestCase):

    def test_initialize(self):
        _ic = incar()
        _ic.parse_tags(encutPw=100)
        _ic.parse_tags(ENCUT=100)

class test_tag_manipulation(ut.TestCase):

    def test_parse_tag_and_tag_vals(self):
        _ic = incar(encutPw=100, ENCUTGW=50, gga="PE", ISTART=1)
        # get the same value by either INCAR tag or its "n a" equivalent
        self.assertListEqual([100], _ic.tag_vals("ENCUT"))
        self.assertListEqual([100], _ic.tag_vals("encutPw"))
        self.assertListEqual([50], _ic.tag_vals("encutPwGw"))
        self.assertListEqual([50], _ic.tag_vals("ENCUTGW"))
        self.assertListEqual([1], _ic.tag_vals("ISTART"))
        self.assertListEqual([1], _ic.tag_vals("restartWave"))
        self.assertListEqual(["PE"], _ic.tag_vals("GGA"))
        self.assertListEqual(["PE"], _ic.tag_vals("gga"))
        _ic.parse_tags(ALGO=True)
        self.assertListEqual([True], _ic.tag_vals("ALGO"))
        _ic.parse_tags(scfAlgo=False)
        self.assertListEqual([False], _ic.tag_vals("scfAlgo"))
        # parse invalid tag
        _ic.parse_tags(abc=1)
        self.assertListEqual([None], _ic.tag_vals("abc"))
        _ic.parse_tags()


    def test_pop_and_delete_tags(self):
        _ic = incar(ENCUT=100, ENCUTGW=50, gga="PE", ISTART=1, AEXX=0.1, EDIFF=1.0E-6)
        ecut = _ic.pop_tags("ENCUT")[0]
        self.assertEqual(100, ecut)
        self.assertListEqual([None], _ic.tag_vals("ENCUT"))
        _ic.delete_tags("gga","ISTART")
        self.assertListEqual([None,]*2, _ic.tag_vals("GGA", "restartWave"))



class test_read_from_file(ut.TestCase):
    pass



def _verify_incar_from_json(tc, ic, pathJson):
    pass


if __name__ == '__main__':
    ut.main()