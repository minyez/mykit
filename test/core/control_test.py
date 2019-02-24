#!/usr/bin/env python3
# coding = utf8

import unittest as ut
from mykit.core._control import control_map


class test_control_map(ut.TestCase):

    def test_tags_mapping(self):

        __mapDict = {
            "abc": {"l":"abc", "a": "123", "U": "ABC"},
            "def": {"l":"def", "a": "456", "U": "DEF"},
            "ghi": {"l":"ghi", "a": "789", "U": "GHI"},
        }
        __mappedTag = control_map._tags_mapping(__mapDict, "l", "a", "abc", "def", "ghi")
        self.assertTupleEqual(("123","456","789"), __mappedTag) 
        __mappedTag = control_map._tags_mapping(__mapDict, "a", "l", "abc", "def", "ghi")
        self.assertTupleEqual((None,)*3, __mappedTag) 
        __mappedTag = control_map._tags_mapping(__mapDict, "a", "l", "123", "456", "789")
        self.assertTupleEqual(("abc","def","ghi"), __mappedTag) 
        __mappedTag = control_map._tags_mapping(__mapDict, "l", "U", "abc", "def", "ghi")
        self.assertTupleEqual(("ABC","DEF","GHI"), __mappedTag) 
        __mappedTag = control_map._tags_mapping(__mapDict, "U", "l", "abc", "def", "ghi")
        self.assertTupleEqual((None,)*3, __mappedTag) 

if __name__ == '__main__':
    ut.main()