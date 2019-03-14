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
        ("A-B-C D-E", ["A", "B", "B", "C", "D", "E"]),
        ("GM-W-X-L-GM-X", ["GM", "W", "W", "X", "X", "L", "L", "GM", "GM", "X"]),
    )
    testSetsBadKlineInvalid = (
        ("GMX-W", "GMX-W" ),
        ("GM-W LDL-L", "LDL-L"),
        ("GM-W LDL-", "LDL-"),
        ("GM-W LD--L", "LD--L"),
    )
    testSetsBadKlineZero = (
        ("GM-GM", "GM-GM" ),
        ("GM-W LD-LD", "LD-LD"),
        ("GM-W-W L-LD", "W-W"),
    )
    testSetsBadKsymInvalid = (
        (["GM", "W", "LDL", "L"], "LDL"),
        (["GMX", "W"], "GMX"),
    )
    testSetsBadKsymOddLen = (
        (["GM", "W", "LD"], 3),
        (["GM", "W", "L", "X", "GM"], 5),
    )


    def test_decoder(self):
        for _i, (kline, ksyms) in enumerate(self.testSetsValid):
            self.assertListEqual(kpath_decoder(kline), ksyms)
        for _i, (kline, badPart) in enumerate(self.testSetsBadKlineInvalid):
            self.assertRaisesRegex(KmeshError, \
                r"Invalid kpath line string: {}".format(badPart),\
                kpath_decoder, kline)
        for _i, (kline, badPart) in enumerate(self.testSetsBadKlineZero):
            self.assertRaisesRegex(KmeshError, \
                r"kpath with zero length: {}".format(badPart),\
                kpath_decoder, kline)

    def test_encoder(self):
        for _i, (kline, ksyms) in enumerate(self.testSetsValid):
            self.assertEqual(kpath_encoder(ksyms), kline)
        for _i, (ksyms, badPart) in enumerate(self.testSetsBadKsymInvalid):
            self.assertRaisesRegex(KmeshError, \
                r"Invalid kpoint symbol: {}".format(badPart),\
                kpath_encoder, ksyms)
        for _i, (ksyms, badPart) in enumerate(self.testSetsBadKsymOddLen):
            self.assertRaisesRegex(KmeshError, \
                r"require even length, received {}".format(badPart),\
                kpath_encoder, ksyms)
        


if __name__ == '__main__':
    ut.main()
