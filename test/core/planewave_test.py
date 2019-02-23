#!/usr/bin/env python3
# coding = utf-8

import unittest as ut
from mykit.core.planewave import plane_wave, planewaveError

# These test tuples should have the same length
haveImplePwTags = ("encutPw", "encutPwGw", "restartWave")
haveImpleVaspMap = ("ENCUT", "ENCUTGW", "ISTART")
haveImpleQeMap = (None, None, None)
haveImpleAbiMap = (None, None, None)

class test_pwtags(ut.TestCase):

    pass

    def test_implemented_pwtags(self):
        implePwTags = plane_wave.map2pwtags(getAll=True)
        for _t in haveImplePwTags:
            self.assertIn(_t, implePwTags)


class test_tag_mapping(ut.TestCase):

    pass

    def test_map_tags(self):
        '''Test the map_tags classmethod
        '''
        for _i,_v in enumerate(haveImplePwTags):
            # single tag
            self.assertTupleEqual((haveImplePwTags[_i],), plane_wave.map_tags_in_pw(haveImplePwTags[_i]))
            self.assertTupleEqual((haveImpleVaspMap[_i],), plane_wave.map_tags_in_pw(haveImplePwTags[_i], progTo="vasp"))
            self.assertTupleEqual((haveImpleQeMap[_i],), plane_wave.map_tags_in_pw(haveImplePwTags[_i], progTo="qe"))
            # multiple tags
            self.assertTupleEqual(haveImplePwTags[:_i+1], plane_wave.map_tags_in_pw(*haveImplePwTags[:_i+1]))
            self.assertTupleEqual(haveImpleVaspMap[:_i+1], plane_wave.map_tags_in_pw(*haveImplePwTags[:_i+1], progTo="vasp"))
            self.assertTupleEqual(haveImpleQeMap[:_i+1], plane_wave.map_tags_in_pw(*haveImplePwTags[:_i+1], progTo="qe"))

        # Map to plane_wave tags from tags for particular program. 
        # Should be noticed that the tag list to map from should not have None in it, 
        # otherwise the following test must fail.
        self.assertTupleEqual(haveImplePwTags, plane_wave.map_tags_in_pw(*haveImpleVaspMap, progFrom="vasp"))
        # self.assertTupleEqual(self.__haveImplePwTags, plane_wave.map_tags(*self.__haveImpleQeMap, progFrom="qe"))

        # Map tags from one program to another.
        # Should be noticed that the tag list to map from should not have None in it, 
        # otherwise the following test must fail.
        self.assertTupleEqual(haveImpleQeMap, plane_wave.map_tags_in_pw(*haveImpleVaspMap, progFrom="vasp", progTo="qe"))
    
    def test_map2pwtags(self):
        '''Test the map2pwtags classmethod
        '''
        self.assertTupleEqual(haveImplePwTags, plane_wave.map2pwtags(*haveImplePwTags, progFrom="n a"))
        
    def test_map_from_pwtags(self):
        '''Test the map_from_pwtags classmethod
        '''
        self.assertTupleEqual(haveImpleVaspMap, plane_wave.map_from_pwtags(*haveImplePwTags, progTo="vasp"))
        

class test_tag_manipulation(ut.TestCase):

    pass
    
    def test_initialization(self):
        # empty initialization
        _pw = plane_wave()
        _pw = plane_wave(encutPw=300, restartWave=0)
        self.assertListEqual([300, 0], _pw.tag_vals("encutPw", "restartWave"))
        self.assertListEqual([300, 0], _pw.tag_vals("ENCUT", "restartWave", progName="vasp"))

    def test_parse_tags(self):
        _pw = plane_wave(restartWave=0)
        self.assertListEqual([0,], _pw.tag_vals("restartWave"))
        _pw.parse_tags(encutPw=200)
        self.assertListEqual([200,], _pw.tag_vals("encutPw"))
        

if __name__ == "__main__":
    # check consistency of input maps
    maps = (haveImplePwTags,haveImpleAbiMap, haveImpleQeMap, haveImpleVaspMap)
    mapItems = tuple(map(len, maps))
    try:
        assert all(i == mapItems[0] for i in mapItems)
    except AssertionError:
        raise IndexError("Test maps are of different length.")

    ut.main()
