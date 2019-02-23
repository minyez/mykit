#!/usr/bin/env python3
# coding = utf-8

import unittest as ut
from mykit.vasp.incar import incar, incarError

class test_direct_set(ut.TestCase):

    def test_initialize(self):
        _ic = incar()
        _ic.parse_tags(encutPw=100)
        self.assertEqual([100], _ic.tag_vals("encutPw"))
        # self.assertEqual([100], _ic.tag_vals("INCAR"))
        _ic.parse_tags(INCAR=100)
        self.assertEqual([100], _ic.tag_vals("encutPw"))
        # self.assertEqual([100], _ic.tag_vals("INCAR"))


class test_read_from_file(ut.TestCase):
    pass



def _verify_incar_from_json(tc, ic, pathJson):
    pass


if __name__ == '__main__':
    ut.main()