#!/usr/bin/env python3
# coding = utf-8

import json
import os
import re
import tempfile
import unittest as ut

from mykit.vasp.incar import (Incar, IncarError, _get_para_tags_from_nproc,
                              _get_xc_tags_from_xcname)


class test_tag_manipulation(ut.TestCase):

    def test_empty_init(self):
        _ic = Incar()
        _ic.parse_tags(encutPw=100)
        _ic.parse_tags(ENCUT=100)

    def test_parse_tag_and_tag_vals(self):
        _ic = Incar(encutPw=100, ENCUTGW=50, gga="PE", ISTART=1)
        # get the same value by either INCAR tag or its "n a" equivalent
        self.assertListEqual([100], _ic.tag_vals("ENCUT"))
        self.assertListEqual([100], _ic.tag_vals("encutPw"))
        # self.assertListEqual([50], _ic.tag_vals("encutPwGw"))
        self.assertListEqual([50], _ic.tag_vals("ENCUTGW"))
        self.assertListEqual([1], _ic.tag_vals("ISTART"))
        # self.assertListEqual([1], _ic.tag_vals("restartWave"))
        self.assertListEqual(["PE",], _ic.tag_vals("GGA"))
        self.assertListEqual(["PE",], _ic.tag_vals("gga"))
        _ic.parse_tags(ALGO=True)
        self.assertListEqual([True], _ic.tag_vals("ALGO"))
        # _ic.parse_tags(scfAlgo=False)
        # self.assertListEqual([False], _ic.tag_vals("scfAlgo"))
        # parse invalid tag
        _ic.parse_tags(abc=1)
        self.assertListEqual([None,], _ic.tag_vals("abc"))
        _ic.parse_tags()
        _ic.parse_tags(IBRION=1)
        self.assertListEqual([1], _ic.tag_vals("IBRION"))
        self.assertListEqual([1], _ic.tag_vals("ionAlgo"))
        _ic.parse_tags(ionAlgo=2)
        self.assertListEqual([2], _ic.tag_vals("IBRION"))
        
        # prefer to INCAR tags to mykit tags, even the mykit tag is parsed later
        _ic.parse_tags(GGA="CA", gga="PE")
        self.assertListEqual(["CA"], _ic.tag_vals("GGA"))
        _ic.parse_tags(IBRION=3, ionAlgo=4)
        self.assertListEqual([3], _ic.tag_vals("IBRION"))

    def test_pop_and_delete_tags(self):
        _ic = Incar(ENCUT=100, ENCUTGW=50, gga="PE", ISTART=1, AEXX=0.1, EDIFF=1.0E-6)
        ecut = _ic.pop_tags("ENCUT")[0]
        self.assertEqual(100, ecut)
        self.assertListEqual([None], _ic.tag_vals("ENCUT"))
        _ic.delete_tags("gga","ISTART")
        self.assertListEqual([None,]*2, _ic.tag_vals("GGA", "restartWave"))

    def test_getitem(self):
        _ic = Incar(ENCUT=100, ENCUTGW=50, gga="PE", ISTART=1, AEXX=0.1, EDIFF=1.0E-6)
        self.assertEqual(100, _ic["ENCUT"])
        self.assertEqual(50, _ic["ENCUTGW"])
        self.assertEqual("PE", _ic["GGA"])
        self.assertEqual(1, _ic["ISTART"])
        self.assertRaises(IncarError, _ic.__getitem__, "gga")
    
    def test_setitem(self):
        _ic = Incar()
        _ic["ENCUT"] = 100
        self.assertListEqual([100,], _ic.tag_vals("ENCUT"))
        _ic["GGA"] = "PE"
        self.assertListEqual(["PE",], _ic.tag_vals("GGA"))
        self.assertRaises(IncarError, _ic.__setitem__, "gga", "PE")


class test_Incar_io(ut.TestCase):

    def test_read_from_file(self):
        _Incar = None
        _countGood = 0
        _countBad = 0
        _countVerified = 0
        __IncarDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'testdata', 'vasp')
        if os.path.isdir(__IncarDir):
            for _f in os.listdir(__IncarDir):
                if re.match('^INCAR_[0-9]+$', _f):
                    _countGood += 1
                    _i = _f.split('_')[1]
                    _path = os.path.join(__IncarDir, _f)
                    _Incar = Incar.read_from_file(_path)
                    _verifyJson = os.path.join(__IncarDir, 'verify_incar_'+_i+'.json')
                    if os.path.isfile(_verifyJson):
                        # print(_Incar)
                        _vs = _verify_Incar_from_json(self, _Incar, _verifyJson)
                        if _vs:
                            _countVerified += 1
                if re.match('^bad_INCAR_[0-9]+$', _f):
                    _countBad += 1
                    _i = _f.split('_')[2]
                    _path = os.path.join(__IncarDir, _f)
                    self.assertRaises(IncarError, Incar.read_from_file, _path)
        print("{} good INCARs readed ({} verified by JSON file). {} bad INCARs raised.".format(_countGood, _countVerified, _countBad))

    def test_print_write(self):
        _ic = Incar(comment="Test", ENCUT=200.0)
        self.assertEqual("Test\nENCUT = 200.0", _ic.__str__())
        tf = tempfile.NamedTemporaryFile()
        _ic.write(pathIncar=tf.name)
        tf.close()


class test_Incar_factory(ut.TestCase):

    def test_minimal_incar(self):
        self.assertRaisesRegex(IncarError, \
            r"Task name not supported: unknown-task. Should be *", \
            Incar.minimal_incar, "unknown-task")
        Incar.minimal_incar("scf")
        Incar.minimal_incar("opt")
        Incar.minimal_incar("dos")
        Incar.minimal_incar("band")
        Incar.minimal_incar("diag")


class test_tag_functions(ut.TestCase):

    def test_get_xc_tags(self):

        xctags = _get_xc_tags_from_xcname("lda")
        self.assertDictEqual(xctags, {"GGA": "CA"})
        xctags = _get_xc_tags_from_xcname("pbe")
        self.assertDictEqual(xctags, {"GGA": "PE"})
        xctags = _get_xc_tags_from_xcname("PBE0")
        self.assertDictEqual(xctags, \
            {"GGA": "PE", "ALGO": "ALL", "LHFCALC": True, "TIME": 0.4, "PRECFOCK": "Normal"})
        xctags = _get_xc_tags_from_xcname("HSE06")
        self.assertDictEqual(xctags, \
            {"GGA": "PE", "ALGO": "ALL", "LHFCALC": True, "TIME": 0.4, "PRECFOCK": "Normal", "HFSCREEN": 0.2})
        xctags = _get_xc_tags_from_xcname("HF")
        self.assertDictEqual(xctags, \
            {"GGA": "PE", "ALGO": "ALL", "LHFCALC": True, "TIME": 0.4, "PRECFOCK": "Normal", "AEXX": 1.0})
        self.assertRaisesRegex(IncarError, r"XC name is not supported: *", _get_xc_tags_from_xcname, "non-xc", ignore_error=False)

    def test_get_para_tags(self):

        paratags = _get_para_tags_from_nproc(1)
        self.assertDictEqual(paratags, {})
        paratags = _get_para_tags_from_nproc(4)
        self.assertDictEqual(paratags, {"NPAR": 2, "KPAR": 2})
        paratags = _get_para_tags_from_nproc(27)
        self.assertDictEqual(paratags, {"NPAR": 9, "KPAR": 3})
        

def _verify_Incar_from_json(tc, ic, pathJson):
    assert isinstance(tc, ut.TestCase)
    assert isinstance(ic, Incar)
    with open(pathJson, 'r') as f:
        _vDict = json.load(f)
    verifyMsg = "INCAR verification failed: {}".format(pathJson)
    tc.assertDictEqual(_vDict, ic.tags, msg=verifyMsg)
    return True


if __name__ == '__main__':
    ut.main()
