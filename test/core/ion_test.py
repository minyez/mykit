#!/usr/bin/env python3
# coding = utf-8

import unittest as ut
from mykit.core.ion import IonError
from mykit.core.ion import ion_control as ionc


# These test tuples should have the same length
haveImpleIonTags = ("nIonStepsMax", "convIon", "ionAlgo")
haveImpleVaspMap = ("NSW", "EDIFFG","IBRION")
haveImpleQeMap = (None,None,None)
haveImpleAbiMap = (None,None,None)

class test_iontags(ut.TestCase):

    def test_implemented_iontags(self):
        impleIonTags = ionc.map_to_mykit_tags(getAll=True)
        for _t in haveImpleIonTags:
            self.assertIn(_t, impleIonTags)

class test_tag_mapping(ut.TestCase):

    pass

    def test_map_tags_in_ion(self):

        for i, _v in enumerate(haveImpleIonTags):
            # single tag
            self.assertTupleEqual((haveImpleIonTags[i],), ionc.map_tags(haveImpleIonTags[i]))
            self.assertTupleEqual((haveImpleVaspMap[i],), ionc.map_tags(haveImpleIonTags[i], progTo="vasp"))
            self.assertTupleEqual((haveImpleQeMap[i],), ionc.map_tags(haveImpleIonTags[i], progTo="qe"))
            # multiple tags
            self.assertTupleEqual(haveImpleIonTags[:i+1], ionc.map_tags(*haveImpleIonTags[:i+1]))
            self.assertTupleEqual(haveImpleVaspMap[:i+1], ionc.map_tags(*haveImpleIonTags[:i+1], progTo="vasp"))
            self.assertTupleEqual(haveImpleQeMap[:i+1], ionc.map_tags(*haveImpleIonTags[:i+1], progTo="qe"))

if __name__ == '__main__':
    ut.main()