#!/usr/bin/env python3
# coding = utf-8

import unittest as ut
from mykit.vasp.incar import incar, incarError

class test_direct_set(ut.TestCase):

    def test_initialize(self):
        _ic = incar()
        _ic.parse_tags(encutPw=100)
        _ic.parse_tags(ENCUT=100)
    
    def test_tag_vals(self):
        # _ic = incar(ENCUT=100, ENCUTGW=50, GGA="PE")
        _ic = incar(encutPw=100, ENCUTGW=50, gga="PE", ISTART=1)
        # get the same value by either INCAR tag or its "n a" equivalent
        self.assertListEqual([100], _ic.tag_vals("ENCUT"))
        self.assertListEqual([100], _ic.tag_vals("encutPw"))
        self.assertListEqual([50], _ic.tag_vals("encutPwGw"))
        self.assertListEqual([50], _ic.tag_vals("ENCUTGW"))
        self.assertListEqual(["PE"], _ic.tag_vals("GGA"))
        self.assertListEqual(["PE"], _ic.tag_vals("gga"))
        # parse invalid tag
        _ic.parse_tags(abc=1)
        self.assertListEqual([None], _ic.tag_vals("abc"))
        _ic.parse_tags()
        


class test_read_from_file(ut.TestCase):
    pass



def _verify_incar_from_json(tc, ic, pathJson):
    pass


if __name__ == '__main__':
    ut.main()