#!/usr/bin/env python3
# coding = utf-8

import unittest as ut

from mykit.vasp.kpoints import KpointsError, kpoints


class test_kpoints_init(ut.TestCase):

    def test_raise_at_initiate(self):
        self.assertRaisesRegex(KpointsError, \
            r"kdense must not be negative", \
            kpoints, kdense=-1)
        self.assertRaisesRegex(KpointsError, \
            r"Enter at least one of *", \
            kpoints)
        self.assertRaisesRegex(KpointsError, \
            r"Unknown KPOINTS mode: *", \
            kpoints, kmode="D", kgrid=[6,6,6])
        self.assertRaisesRegex(KpointsError, \
            r"kpath should be specified for line mode", \
            kpoints, kmode="L", kdense=10)
        self.assertRaisesRegex(KpointsError, \
            r"Fully automatic mode needs positive kdense*", \
            kpoints, kmode="A", kgrid=[6,6,6])
        self.assertRaisesRegex(KpointsError, \
            r"kgrid should be specified*", \
            kpoints, kmode="G", kdense=10)
        self.assertRaisesRegex(KpointsError, \
            r"kgrid should be specified*", \
            kpoints, kmode="M", kdense=10)

    def test_auto_mode_check(self):
        # use G-mode when kgrid is specified
        kp = kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kgrid=[5,5,5])
        self.assertEqual("G", kp.mode)
        # when kdense is positive and kmode is not "L", switch to automatic
        kp = kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kdense=10)
        self.assertEqual("A", kp.mode)
        # when kpath is specified, always use line mode
        _kpath = [(0.0,0.0,0.0), (0.0,0.5,0.0)]
        kp = kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kdense=10, kpath=_kpath)
        self.assertEqual("L", kp.mode)

    def test_raise_in_G_M_mode(self):
        kp = kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kgrid=[5,5], kmode="G")
        self.assertRaisesRegex(KpointsError, \
            r"Bad kgrid format for G/M mode: *", \
            print, kp)
    
    def test_raise_in_L_mode(self):
        _kpathsBadDiffLen = [
                {
                    "coordinates" :[(0.0, 0.0, 0.0)],
                    "symbols": ["GM", "L"],
                },
                {
                    "coordinates" :[(0.0, 0.0, 0.0),(0.5,0.5,0.5)],
                    "symbols": ["GM",],
                },
            ]
        _kpathsBadOdd = [
                {
                    "coordinates" :[(0.0, 0.0, 0.0)],
                    "symbols": ["GM"],
                },
                {
                    "coordinates" :[(0.0, 0.0, 0.0),(0.5,0.5,0.5),(0.5,0.5,0.5)],
                    "symbols": ["GM", "L", "L"],
                },
            ]
        for _kpath in _kpathsBadDiffLen:
            kp = kpoints(comment="KPOINTS for test_raise_in_L_mode", \
                kmode="L", kdense=15, kpath=_kpath)
            self.assertRaisesRegex(KpointsError, \
                "Inconsistent length of symbols and coordiantes", \
                print, kp)
        for _kpath in _kpathsBadOdd:
            kp = kpoints(comment="KPOINTS for test_raise_in_L_mode", \
                kmode="L", kdense=15, kpath=_kpath)
            self.assertRaisesRegex(KpointsError, \
                "Odd length found for symbols/coordiantes, require even.", \
                print, kp)


if __name__ == '__main__':
    ut.main()
