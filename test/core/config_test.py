#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import sys
import tempfile
import unittest as ut
from shutil import copy2, move, which

from mykit.core.config import ConfigError, global_config


class check_glocal_config(ut.TestCase):

    def test_non_existent_option(self):
        
        os.unsetenv(global_config.env_var())
        _c = global_config()
        self.assertEqual(None, _c.get())
        self.assertEqual('', global_config.get_doc('opt_not_exist1'))
        self.assertRaisesRegex(ConfigError, r"Unknown option: opt_not_exist1", _c.get, 'opt_not_exist1')
        self.assertRaisesRegex(ConfigError, r"Unknown option: *", _c.get, 'opt_not_exist1', 'opt_not_exist2')

    def test_defaults(self):
        # back up the custom JSON if there is
        _path = global_config._get_dejson_path()
        _hasCustom = os.path.isfile(_path)
        if _hasCustom:
            _tf = tempfile.NamedTemporaryFile()
            move(_path, _tf.name)

        _c = global_config()
        self.assertEqual(which('mpirun'), _c.get('mpiExec'))
        self.assertEqual((which('vasp_std'), which('mpirun')), \
            _c.get('vaspStdExec','mpiExec'))
        self.assertEqual(_c.get('mpiExec', doc=True), "the MPI executable to use")
        self.assertEqual(_c.get('mpiExec', 'vaspStdExec', doc=True), \
            ("the MPI executable to use", "Path of `vasp_std` executive"))

        if _hasCustom:
            copy2(_tf.name, _path)
            _tf.close()

    def test_load_custom(self):
        '''test custom configuration file created temporarily'''
        _opts = {"vaspStdExec": "vasp", "vaspPawPbe": "/Users/pbe", "vaspPawLda": "/Users/lda"}

        _tf = tempfile.NamedTemporaryFile()
        with open(_tf.name, 'w') as _f:
            json.dump(_opts, _f)
        # print(_tf.name)
        os.environ[global_config.env_var()] = _tf.name

        _c = global_config()
        # os.unsetenv(global_config.env_var())
        del os.environ[global_config.env_var()]
        # print(os.environ[global_config.env_var()])
        # global_config.print_opts()
        self.assertEqual(_c.get('vaspStdExec', doc=True), "Path of `vasp_std` executive")
        self.assertEqual("vasp", _c.get('vaspStdExec'))
        self.assertTupleEqual(("/Users/pbe", "/Users/lda"), _c.get('vaspPawPbe', 'vaspPawLda'))
        # doc string should not be touch
        _tf.close()

    # def test_print(self):
    #     global_config.print_opts

if __name__ == "__main__":
    ut.main()
