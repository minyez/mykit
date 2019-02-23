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

    def test_implemented_pwtags(self):
        implePwTags = plane_wave.map2PwTags()
        for _t in haveImplePwTags:
            self.assertIn(_t, implePwTags)

class test_tag_mapping(ut.TestCase):

    def test_map_tags(self):
        '''Test the map_tags classmethod
        '''
        self.assertTupleEqual(haveImplePwTags, plane_wave.map_tags(*haveImplePwTags))
        # Map plane_wave tags to tags for particular program. 
        self.assertTupleEqual(haveImpleQeMap, plane_wave.map_tags(*haveImplePwTags, progTo="qe"))
        self.assertTupleEqual(haveImpleVaspMap, plane_wave.map_tags(*haveImplePwTags, progTo="vasp"))

        # Map to plane_wave tags from tags for particular program. 
        # Should be noticed that the tag list to map from should not have None in it, 
        # otherwise the following test must fail.
        self.assertTupleEqual(haveImplePwTags, plane_wave.map_tags(*haveImpleVaspMap, progFrom="vasp"))
        # self.assertTupleEqual(self.__haveImplePwTags, plane_wave.map_tags(*self.__haveImpleQeMap, progFrom="qe"))

        # Map tags from one program to another.
        # Should be noticed that the tag list to map from should not have None in it, 
        # otherwise the following test must fail.
        self.assertTupleEqual(haveImpleQeMap, plane_wave.map_tags(*haveImpleVaspMap, progFrom="vasp", progTo="qe"))

class test_intialize(ut.TestCase):
    
    def test_direct_init(self):
        _pw = plane_wave(encutPw=300, restartWave=0)
        self.assertListEqual([300, 0], _pw.tag_vals("encutPw", "restartWave"))
        

if __name__ == "__main__":
    # check consistency of input maps
    maps = (haveImplePwTags,haveImpleAbiMap, haveImpleQeMap, haveImpleVaspMap)
    mapItems = tuple(map(len, maps))
    try:
        assert all(i == mapItems[0] for i in mapItems)
    except AssertionError:
        raise IndexError("Test maps are of different length.")

    ut.main()
