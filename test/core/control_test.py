#!/usr/bin/env python3
# coding = utf8

import unittest as ut
from mykit.core._control import tags_mapping, parse_to_tagdict, check_valid_map


class test_control_map(ut.TestCase):

    _mapDict = {
        "abc": {"l":"abc", "C": "Abc", "U": "ABC"},
        "def": {"l":"def", "C": "Def", "U": "DEF"},
        "ghi": {"l":"ghi", "C": "Ghi", "U": "GHI"},
    }

    def test_tags_mapping(self):
        _mappedTag = tags_mapping(self._mapDict, "l", "C", "abc", "def", "ghi")
        self.assertTupleEqual(("Abc","Def","Ghi"), _mappedTag) 
        _mappedTag = tags_mapping(self._mapDict, "C", "l", "abc", "def", "ghi")
        self.assertTupleEqual((None,)*3, _mappedTag) 
        _mappedTag = tags_mapping(self._mapDict, "C", "l", "Abc", "Def", "Ghi")
        self.assertTupleEqual(("abc","def","ghi"), _mappedTag) 
        _mappedTag = tags_mapping(self._mapDict, "l", "U", "abc", "def", "ghi")
        self.assertTupleEqual(("ABC","DEF","GHI"), _mappedTag) 
        _mappedTag = tags_mapping(self._mapDict, "U", "l", "abc", "def", "ghi")
        self.assertTupleEqual((None,)*3, _mappedTag) 
    
    def test_check_valid_map(self):
        self.assertTrue({"a":{}, "b":{}})
        self.assertTrue(self._mapDict, "l")
    
    def test_parse_to_tagdict(self):
        tags = {}
        parse_to_tagdict(tags, self._mapDict, "l", abc=1)
        self.assertDictEqual({"abc":1}, tags)
        parse_to_tagdict(tags, self._mapDict, "C", Abc=3)
        self.assertDictEqual({"abc":3}, tags)
        # No change when parsing an unknown tag
        parse_to_tagdict(tags, self._mapDict, "l", noexist=1)
        self.assertDictEqual({"abc":3}, tags)
        # No change when parsing a tag that does not correspond to its correct progName
        parse_to_tagdict(tags, self._mapDict, "l", Abc=2)
        self.assertDictEqual({"abc":3}, tags)
        parse_to_tagdict(tags, self._mapDict, "U", ABC=2)
        self.assertDictEqual({"abc":2}, tags)
        # parse multiple tags
        parse_to_tagdict(tags, self._mapDict, "U", DEF=1, GHI=3)
        self.assertDictEqual({"abc":2, "def":1, "ghi":3}, tags)

    def test_extract_from_tagdict(self):
        pass
        

if __name__ == '__main__':
    ut.main()