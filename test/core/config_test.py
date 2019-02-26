#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest as ut
import tempfile
import json
import os
import sys
from shutil import copy2, which, move
from mykit.core.config import global_config

class check_glocal_config(ut.TestCase):

    def test_non_existent_option(self):
        self.assertEqual(None, global_config.get())
        self.assertEqual('', global_config.get(doc=True))
        self.assertEqual(None, global_config.get('opt_not_exist1'))
        self.assertEqual((None, None), global_config.get('opt_not_exist1', 'opt_not_exist2'))

    def test_defaults(self):
        os.unsetenv(global_config._env_var())
        # back up the custom JSON if there is
        _path = global_config._get_dejson_path()
        _hasCustom = os.path.isfile(_path)
        if _hasCustom:
            _tf = tempfile.NamedTemporaryFile()
            move(_path, _tf.name)

        self.assertEqual(which('mpirun'), global_config.get('mpiExec'))
        self.assertEqual((which('vasp_std'), which('mpirun')), \
            global_config.get('vaspStdExec','mpiExec'))
        self.assertEqual(global_config.get('mpiExec', doc=True), "the MPI executable to use")
        self.assertEqual(global_config.get('mpiExec', 'vaspStdExec', doc=True), \
            ("the MPI executable to use", "Path of `vasp_std` executive"))

        if _hasCustom:
            copy2(_tf.name, _path)
            _tf.close()

    def test_load_custom(self):
        '''test custom configuration file created temporarily'''
        _envVar = global_config._env_var()
        _opts = {"vaspStdExec": "vasp"}

        _tf = tempfile.NamedTemporaryFile()
        with open(_tf.name, 'w') as _f:
            json.dump(_opts, _f, indent=2)

        os.environ[_envVar] = _tf.name
        self.assertEqual("vasp", global_config.get('vaspStdExec'))
        self.assertEqual(global_config.get('vaspStdExec', doc=True), "Path of `vasp_std` executive")
        # doc string should not be touch
        _tf.close()

    def test_print(self):
        global_config.print_opts

if __name__ == "__main__":
    ut.main()