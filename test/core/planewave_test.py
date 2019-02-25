#!/usr/bin/env python3
# coding = utf-8

import unittest as ut
from mykit.core.planewave import planewave_control as pwc
from mykit.core.planewave import planewaveError

# These test tuples should have the same length
haveImplePwTags = ("encutPw", "encutPwGw", "restartWave")
haveImpleVaspMap = ("ENCUT", "ENCUTGW", "ISTART")
haveImpleQeMap = (None, None, None)
haveImpleAbiMap = (None, None, None)

class test_pwtags(ut.TestCase):

    pass

    def test_implemented_pwtags(self):
        implePwTags = pwc.map2pwtags(getAll=True)
        for _t in haveImplePwTags:
            self.assertIn(_t, implePwTags)


class test_tag_mapping(ut.TestCase):

    pass

    def test_map_tags_in_pw(self):
        '''Test the map_tags classmethod
        '''
        for _i,_v in enumerate(haveImplePwTags):
            # single tag
            self.assertTupleEqual((haveImplePwTags[_i],), pwc.map_tags_in_pw(haveImplePwTags[_i]))
            self.assertTupleEqual((haveImpleVaspMap[_i],),pwc.map_tags_in_pw(haveImplePwTags[_i], progTo="vasp"))
            self.assertTupleEqual((haveImpleQeMap[_i],), pwc.map_tags_in_pw(haveImplePwTags[_i], progTo="qe"))
            # multiple tags
            self.assertTupleEqual(haveImplePwTags[:_i+1], pwc.map_tags_in_pw(*haveImplePwTags[:_i+1]))
            self.assertTupleEqual(haveImpleVaspMap[:_i+1],pwc.map_tags_in_pw(*haveImplePwTags[:_i+1], progTo="vasp"))
            self.assertTupleEqual(haveImpleQeMap[:_i+1], pwc.map_tags_in_pw(*haveImplePwTags[:_i+1], progTo="qe"))

        # Map to plane_wave tags from tags for particular program. 
        # Should be noticed that the tag list to map from should not have None in it, 
        # otherwise the following test must fail.
        self.assertTupleEqual(haveImplePwTags, pwc.map_tags_in_pw(*haveImpleVaspMap, progFrom="vasp"))

        # Map tags from one program to another.
        # Should be noticed that the tag list to map from should not have None in it, 
        # otherwise the following test must fail.
        self.assertTupleEqual(haveImpleQeMap, pwc.map_tags_in_pw(*haveImpleVaspMap, progFrom="vasp", progTo="qe"))
    
    def test_map2pwtags(self):
        '''Test the map2pwtags classmethod
        '''
        self.assertTupleEqual(haveImplePwTags, pwc.map2pwtags(*haveImplePwTags, progFrom="n a"))
        
    def test_map_from_pwtags(self):
        '''Test the map_from_pwtags classmethod
        '''
        self.assertTupleEqual(haveImpleVaspMap, pwc.map_from_pwtags(*haveImplePwTags, progTo="vasp"))
        

class test_tag_manipulation(ut.TestCase):

    pass
    
    def test_initialization(self):
        # empty initialization
        _pw = pwc("n a")
        _pw = pwc("n a", encutPw=300, restartWave=0)
        self.assertListEqual([300, 0], _pw.tag_vals("n a", "encutPw", "restartWave"))
        self.assertListEqual([300, 0], _pw.tag_vals("vasp", "ENCUT", "restartWave"))
        self.assertListEqual([300, 0], _pw.tag_vals("vasp", "ENCUT", "ISTART"))

    def test_parse_tags(self):
        _pw = pwc("n a", restartWave=0)
        self.assertListEqual([0], _pw.tag_vals("n a", "restartWave"))
        _pw.parse_tags("n a", encutPw=200)
        self.assertListEqual([200], _pw.tag_vals("n a", "encutPw"))
    
    def test_pop_delete_tags(self):
        _pw = pwc("n a", encutPw=300, restartWave=1)
        self.assertListEqual([300, 1], _pw.pop_tags("n a", "encutPw", "restartWave"))
        self.assertListEqual([None, None], _pw.tag_vals("n a", "encutPw", "restartWave"))
        _pw.parse_tags("n a", encutPw=100)
        self.assertListEqual([100], _pw.tag_vals("vasp", "ENCUT"))
        _pw.delete_tags("vasp", "ENCUT")
        self.assertListEqual([None], _pw.tag_vals("n a", "encutPw"))
        

if __name__ == "__main__":
    # check consistency of input maps
    maps = (haveImplePwTags,haveImpleAbiMap, haveImpleQeMap, haveImpleVaspMap)
    mapItems = tuple(map(len, maps))
    try:
        assert all(i == mapItems[0] for i in mapItems)
    except AssertionError:
        raise IndexError("Test maps are of different length.")

    ut.main()
