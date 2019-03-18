#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tempfile
import unittest as ut

from mykit.core.config import global_config
from mykit.core.log import Verbose


class direct_Verbose_test(ut.TestCase):

    _vb = Verbose()

    def test_property(self):
        self._vb.warnLevel = 1
        self._vb.logLevel = 2
        self.assertTrue(self._vb.warnLevel, 1)
        self.assertTrue(self._vb.logLevel, 2)

    def test_indent_level(self):
        _logPref = self._vb._prefix["log"]
        _warnPref = self._vb._prefix["warn"]
        _indent = self._vb._indent
        self._vb.warnLevel = 0
        self._vb.logLevel = 0

        _tf = tempfile.NamedTemporaryFile()
        with open(_tf.name, 'w') as f:
            self._vb.print_log("Hello", "World", file=f)
            self._vb.print_log("Hello", "World", depth=1, file=f)
            self._vb.print_log("Hello", "World", 1, 2, depth=2, file=f)
            # the following log is not printed, as the level 1 > 0 which is set above
            self._vb.print_log("Hello", "World", level=1, file=f)
            self._vb.print_warn("Hello", "World", file=f)
        with open(_tf.name, 'r') as f:
            _str = f.readline().strip("\n")
            self.assertEqual(_str, ' '.join([_logPref, 'Hello World']))
            _str = f.readline().strip("\n")
            self.assertEqual(_str, ' '.join([_indent * 1 + _logPref, 'Hello World']))
            _str = f.readline().strip("\n")
            self.assertEqual(_str, ' '.join([_indent * 2 + _logPref, 'Hello World', "1", "2"]))
            _str = f.readline().strip("\n")
            self.assertNotEqual(_str, ' '.join([_logPref, 'Hello World']))
            self.assertEqual(_str, ' '.join([_warnPref, 'Hello World']))
        _tf.close()

if __name__ == "__main__":
    ut.main()
