#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
from mykit.core.utils import get_dirpath

class file_and_path(unittest.TestCase):

    def test_get_dirpath(self):
        self.assertEqual("/usr/bin", get_dirpath("/usr/bin"))
        self.assertEqual("/bin", get_dirpath("/bin/ls"))


if __name__ == "__main__":
    unittest.main()