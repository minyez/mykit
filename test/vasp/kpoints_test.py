#!/usr/bin/env python3
# coding = utf-8

import tempfile
import unittest as ut

from mykit.vasp.kpoints import Kpoints, KpointsError


class test_kpoints_init(ut.TestCase):

    def test_raise_at_initiate(self):
        self.assertRaisesRegex(KpointsError, \
            r"kdense must not be negative", \
            Kpoints, kdense=-1)
        self.assertRaisesRegex(KpointsError, \
            r"Enter at least one of *", \
            Kpoints)
        self.assertRaisesRegex(KpointsError, \
            r"Unknown KPOINTS mode: *", \
            Kpoints, kmode="D", kdiv=[6,6,6])
        self.assertRaisesRegex(KpointsError, \
            r"kpath should be specified for line mode", \
            Kpoints, kmode="L", kdense=10)
        self.assertRaisesRegex(KpointsError, \
            r"Fully automatic mode needs positive kdense*", \
            Kpoints, kmode="A", kdiv=[6,6,6])
        self.assertRaisesRegex(KpointsError, \
            r"kdiv should be specified*", \
            Kpoints, kmode="G", kdense=10)
        self.assertRaisesRegex(KpointsError, \
            r"kdiv should be specified*", \
            Kpoints, kmode="M", kdense=10)

    def test_auto_mode_check(self):
        '''Check mode identification. Test print as well
        '''
        tf = tempfile.NamedTemporaryFile()
        # use G-mode when kdiv is specified
        kp = Kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kdiv=[5,5,5])
        self.assertEqual("G", kp.mode)
        kp.write(tf.name)
        # when kdense is positive and kmode is not "L", switch to automatic
        kp = Kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kdense=10)
        self.assertEqual("A", kp.mode)
        kp.write(tf.name)
        # when kpath is specified, always use line mode
        _kpath = [(0.0,0.0,0.0), (0.0,0.5,0.0)]
        self.assertRaisesRegex(TypeError, "kpath must be dictionary.", \
            Kpoints, kpath=_kpath, kdense=10)
        _kpath = {
                "symbols": ["GM", "X"],
                "coordinates": [(0.0,0.0,0.0), (0.0,0.5,0.0)],
                }
        kp = Kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kdense=10, kpath=_kpath)
        self.assertEqual("L", kp.mode)
        kp.write(tf.name)
        tf.close()

    def test_raise_in_G_M_mode(self):
        kp = Kpoints(comment="KPOINTS for test_raise_in_G_M_mode", \
            kdiv=[5,5], kmode="G")
        self.assertRaisesRegex(KpointsError, \
            r"Bad kdiv format for G/M mode: *", \
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
            kp = Kpoints(comment="KPOINTS for test_raise_in_L_mode", \
                kmode="L", kdense=15, kpath=_kpath)
            self.assertRaisesRegex(KpointsError, \
                "Inconsistent length of symbols and coordiantes", \
                print, kp)
        for _kpath in _kpathsBadOdd:
            kp = Kpoints(comment="KPOINTS for test_raise_in_L_mode", \
                kmode="L", kdense=15, kpath=_kpath)
            self.assertRaisesRegex(KpointsError, \
                "Odd length found for symbols/coordiantes, require even.", \
                print, kp)


if __name__ == '__main__':
    ut.main()
