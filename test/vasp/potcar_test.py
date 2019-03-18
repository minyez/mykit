#!/usr/bin/env python3
# coding = utf-8

import json
import os
import tempfile
import unittest as ut
from shutil import copy2, move

from mykit.vasp.potcar import Potcar, PotcarError, PotcarSearch


class test_PotcarSearch(ut.TestCase):


    def test_check_names(self):
        # raise first for no element input
        self.assertRaisesRegex(PotcarError, r"should have at least one element", PotcarSearch)
        pts = PotcarSearch("C")
        self.assertTupleEqual(("C",), pts.names)
        pts = PotcarSearch("C", "H", "O")
        self.assertTupleEqual(("C", "H", "O"), pts.names)
        pts = PotcarSearch("C", usegw=True)
        self.assertTupleEqual(("C_GW",), pts.names)
        pts = PotcarSearch("C_GW", "H", usegw=True)
        self.assertTupleEqual(("C_GW", "H_GW"), pts.names)

    def test_xcdir_not_found_when_export(self):
        from mykit.core.config import global_config

        try:
            del os.environ[global_config.env_var()]
        except KeyError:
            pass
        # back up default json
        _path = global_config._get_dejson_path()
        _hasCustom = os.path.isfile(_path)
        if _hasCustom:
            _tf = tempfile.NamedTemporaryFile()
            move(_path, _tf.name)

        _opts = {"vaspPawPbe": "/Users/pawpbe", "vaspPawLda": os.environ["HOME"]}
        with open(_path, 'w') as _f:
            json.dump(_opts, _f, indent=2)
        
        pts = PotcarSearch("C")
        self.assertDictEqual(pts._homePaw, {"PBE": "/Users/pawpbe", "LDA": os.environ["HOME"]})

        try:
            # raise for unknown XC name
            self.assertRaisesRegex(PotcarError, r"XC type not supported: *", \
                pts.export, xc="PW91")
            # raise for paw directory not found
            self.assertRaisesRegex(PotcarError, r"PBE PAW home directory does not exist: *", \
                pts.export, xc="pbe")
            # raise for invalid PAW directory structure
            self.assertRaisesRegex(PotcarError, \
                "POTCAR not found: {}".format(os.path.join(os.environ["HOME"], "C", "POTCAR")), \
                    pts.export, xc="lda")
        except AssertionError as _err:
            # Remove file and copy back the default json
            os.remove(_path)
            if _hasCustom:
                copy2(_tf.name, _path)
                _tf.close()
            raise _err
        else:
            os.remove(_path)
            if _hasCustom:
                copy2(_tf.name, _path)
                _tf.close()

    def test_change_names(self):
        pts = PotcarSearch("C")
        self.assertTupleEqual(pts.names, ("C",))
        pts.names = ("P", "C_GW")
        self.assertTupleEqual(("P", "C_GW"), pts.names)


class test_potcar(ut.TestCase):

    def test_get_enmin_enmax(self):
        '''Test get ENMIN and ENMAX from a fake POTCAR
        '''
        self.assertRaisesRegex(PotcarError, r"POTCAR not found: *", \
            Potcar.get_enmin_enmax, "not-exist-path.POTCAR")
        fakePotPath = os.path.join(os.path.dirname(__file__), \
            '..', 'testdata', 'vasp', 'POTCAR_fake')
        self.assertTupleEqual((310.494, 413.992), Potcar.get_enmin_enmax(fakePotPath))
        badPotPath = os.path.join(os.path.dirname(__file__), \
            '..', 'testdata', 'vasp', 'POTCAR_bad_enm')
        self.assertRaisesRegex(PotcarError, r"Bad POTCAR for ENMAX and ENMIN: *", \
            Potcar.get_enmin_enmax, badPotPath)


if __name__ == '__main__':
    ut.main()
