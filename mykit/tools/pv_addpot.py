#!/usr/bin/env python3
'''Search or create POTCAR file.

The elements can be parsed manually, or read from POSCAR.
'''

import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from mykit.core.utils import trim_after
from mykit.vasp.poscar import poscar
from mykit.vasp.potcar import potcar_search


def pv_addpot():

    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    groupEle = parser.add_mutually_exclusive_group()
    groupEle.add_argument('-e', dest='elements', default=None, nargs='+', \
            help="POTCAR names, e.g. Cu, Fe_sv, H_GW, etc")
    groupEle.add_argument('-i', dest='posin', default='POSCAR', \
            help="POSCAR file to extract elements")
    parser.add_argument("-s", dest='search', action='store_true', \
            help="flag for search mode, searching available POTCARs")
    parser.add_argument("-x", dest='xc', default='PBE', \
            help="XC functional to generate PAW. LDA or PBE (default).")
    parser.add_argument("--gw", dest="usegw", default=False, action='store_true', \
            help="flag for always using GW potential")
    # parser.add_argument"-D", dest='debug', action='store_true', \
    #         help="flag for debug mode")
    opts = parser.parse_args()

    if opts.elements != None:
        ele = opts.elements
    else:
        _pc = poscar.read_from_file(opts.posin)
        ele = _pc.atomTypes
        # Trim ele to get rid of number
        ele = [trim_after(e, r"\d") for e in ele]
    
    pts = potcar_search(*ele, usegw=opts.usegw)
    # if opts.debug:
    #     print(opts.xc)
    #     print(ele)
    if opts.search:
        _home, _dict = pts.search(opts.xc)
        print("Searching POTCARs from: {}".format(_home))
        for name, v in _dict.items():
            print("- {}".format(name))
            for pot, (emin, emax) in v.items():
                print("  {:10s}: {:8.3f} ~{:8.3f}".format(pot, emin, emax))
    else:
        print("Exporting ", *pts.names)
        _home = pts.export(opts.xc)
        print("POTCAR exported from {}".format(_home))
        

if __name__ == "__main__":
    pv_addpot()
