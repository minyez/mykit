#!/usr/bin/env python3
# coding = utf-8

import unittest as ut

from mykit.core.kmesh import (KmeshError, kmesh_control, kpath_decoder,
                              kpath_encoder)


class test_kpath_coder_decoder(ut.TestCase):
    '''Test the behavior of kpath decoder/encoder

    Note:
        The validity of each test set may vary, since `KSYM_PATTERN` in `kmesh`
        might change.
    '''
    testSetsValid = (
        ("A-B-C D-E", [("A", "B"), ("B", "C"), ("D", "E")]),
        ("GM-W-X-L-GM-X", [("GM", "W"), ("W", "X"), ("X", "L"), ("L", "GM"), ("GM", "X")]),
    )
    testSetsBadKline = (
        ("GMX-W", "GMX-W" ),
        ("GM-W LDL-L", "LDL-L"),
        ("GM-W LDL-", "LDL-"),
        ("GM-W LD--L", "LD--L"),
    )
    testSetsBadKsegsSym = (
        ([("GM", "W"), ("LDL", "L")], "LDL"),
        ([("GMX", "W"),], "GMX"),
    )

    def test_decoder(self):
        for _i, (kline, ksegs) in enumerate(self.testSetsValid):
            self.assertListEqual(kpath_decoder(kline), ksegs)
        for _i, (kline, badPart) in enumerate(self.testSetsBadKline):
            self.assertRaisesRegex(KmeshError, \
                r"Invalid kpath line string: {}".format(badPart),\
                kpath_decoder, kline)

    def test_encoder(self):
        for _i, (kline, ksegs) in enumerate(self.testSetsValid):
            self.assertEqual(kpath_encoder(ksegs), kline)
        for _i, (ksegs, badPart) in enumerate(self.testSetsBadKsegsSym):
            self.assertRaisesRegex(KmeshError, \
                r"Invalid kpoint symbol: {}".format(badPart),\
                kpath_encoder, ksegs)


if __name__ == '__main__':
    ut.main()
