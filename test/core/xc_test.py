#!/usr/bin/env python3
# coding = utf8

import unittest as ut
from mykit.core.xc import xc_control as xcc
from mykit.core.xc import xcError

haveImpleXcTags = ("gga", "metagga")
haveImpleVaspMap = ("GGA", "METAGGA")
haveImpleQeMap = (None, None)
haveImpleAbiMap = (None, None)

class test_init(ut.TestCase):

    pass

    def test_direct_init(self):
        _xc = xcc("n a")
        _xc = xcc("n a", gga="PE")
        _xc = xcc("n a", encutGw=500)


class test_map_tags_in_xc(ut.TestCase):
    def test_map_tags_in_xc(self):
        for _i,_v in enumerate(haveImpleXcTags):
            pass
            # single tag
            self.assertTupleEqual((haveImpleXcTags[_i],), xcc.map_tags_in_xc(haveImpleXcTags[_i]))
            self.assertTupleEqual((haveImpleVaspMap[_i],),xcc.map_tags_in_xc(haveImpleXcTags[_i], progTo="vasp"))
            self.assertTupleEqual((haveImpleQeMap[_i],), xcc.map_tags_in_xc(haveImpleXcTags[_i], progTo="qe"))
            # multiple tags
            self.assertTupleEqual(haveImpleXcTags[:_i+1], xcc.map_tags_in_xc(*haveImpleXcTags[:_i+1]))
            self.assertTupleEqual(haveImpleVaspMap[:_i+1],xcc.map_tags_in_xc(*haveImpleXcTags[:_i+1], progTo="vasp"))
            self.assertTupleEqual(haveImpleQeMap[:_i+1], xcc.map_tags_in_xc(*haveImpleXcTags[:_i+1], progTo="qe"))

    def test_map2xctags(self):
        '''Test the map2pwtags classmethod
        '''
        self.assertTupleEqual(haveImpleXcTags, xcc.map2xctags(*haveImpleXcTags, progFrom="n a"))
        
    def test_map_from_xctags(self):
        '''Test the map_from_pwtags classmethod
        '''
        self.assertTupleEqual(haveImpleVaspMap, xcc.map_from_xctags(*haveImpleXcTags, progTo="vasp"))

class test_tag_val(ut.TestCase):

    def test_direct(self):
        _xc = xcc("n a", gga="PE")
        _xc.parse_tags("n a", metagga="SCAN")
        self.assertListEqual(["PE"], _xc.tag_vals("gga"))
        self.assertListEqual(["PE","SCAN"], _xc.tag_vals("gga", "metagga"))


if __name__ == '__main__':
    ut.main()