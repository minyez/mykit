#!/usr/bin/env python3
# coding = utf-8

import unittest as ut

from mykit.vasp.kpoints import KpointsError, kpoints


class test_initiate_from_statement(ut.TestCase):

    def test_raise(self):
        self.assertRaisesRegex(KpointsError, r"Invalid kpoint *", \
            kpoints, kdense=-1)
        self.assertRaisesRegex(KpointsError, r"Enter at least one of *", \
            kpoints)
        self.assertRaisesRegex(KpointsError, r"Fully automatic mode needs positive kdense*", \
            kpoints, kmode="A", kgrid=[6,6,6])
        self.assertRaisesRegex(KpointsError, r"kgrid should be specified*", \
            kpoints, kmode="G", kdense=10)
        self.assertRaisesRegex(KpointsError, r"kgrid should be specified*", \
            kpoints, kmode="M", kdense=10)
        self.assertRaisesRegex(KpointsError, r"Unknown KPOINTS mode:*", \
            kpoints, kmode="D", kgrid=[6,6,6])


if __name__ == '__main__':
    ut.main()
