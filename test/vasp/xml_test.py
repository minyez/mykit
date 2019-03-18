#!/usr/bin/env python3
# coding = utf-8

import os
import unittest as ut

from mykit.vasp.xml import Vasprunxml, VasprunxmlError


class test_vasprunxml_read(ut.TestCase):

    staticXmlId = (0, )

    def test_static_xml(self):
        for id in self.staticXmlId:
            path = os.path.join(os.path.dirname(__file__), '..', \
                'testdata', 'vasprun_{}.xml'.format(id))
            _vxml = Vasprunxml(path)

if __name__ == '__main__':
    ut.main()
