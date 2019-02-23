#!/usr/bin/env python3
# coding = utf-8

import unittest as ut
from mykit.core.planewave import plane_wave

class test_pwtags(ut.TestCase):

    def test_implemented_pwtags(self):
        self.assertTrue(True)
        shouldImplePwTags = ["encutPw", "encutPwGw", "restart"]
        implePwTags = plane_wave.map2PwTags()
        for _t in shouldImplePwTags:
            self.assertIn(_t, implePwTags)


if __name__ == "__main__":
    ut.main()
