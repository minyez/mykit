#!/usr/bin/env python
#
# a simple script to get fundamental gap and k-point-averaged band gap

from __future__ import print_function
from pv_anal_utils import vasp_anal_get_gap,vasp_anal_read_eigen,\
                          vasp_anal_get_kavgap
from argparse import ArgumentParser
import os,sys

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Give the fundamental band gap and k-averaged gap.
    Also these two values between any occupied and empty band is possible
    by adding option: -v vb -c cb
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')
    parser.add_argument("-k",dest='kav',help="flag for print( k-averaged gap",action='store_true')
    parser.add_argument("-v",dest='vb',help="the index of valence band",type=int,default=0)
    parser.add_argument("-c",dest='cb',help="the index of conduction band",type=int,default=0)
    parser.add_argument("-f",dest='fix_k',help="flag for k-point fixing mode for k-averaged gap. 0 for fixing VBM, 1 for fixing CBM ",type=int,default=-1)
    parser.add_argument("-i",dest='flag_inv',help="flag for the inverse of gap.",action='store_true')

    opts = parser.parse_args()
    vb = opts.vb
    cb = opts.cb

# =====================================================

    band_struct = vasp_anal_read_eigen(spinpolarzied=False,debug=opts.debug)
    nelec   = band_struct[0][0]
    nkp     = band_struct[0][1]
    bandmax = band_struct[0][2]

    try:
        assert vb >= 0 and cb >= 0
    except AssertionError:
        print( "Invalid band index")
        sys.exit(1)

    if vb == 0 and cb == 0:
        vb = nelec/2
        cb = nelec/2 + 1
    elif vb == 0:
        if cb <= nelec/2:
            vb = cb
            cb = nelec/2 + 1
        else:
            vb = nelec/2
    elif cb == 0:
        if vb >= nelec/2+1:
            cb = vb
            vb = nelec/2
        else:
            cb = nelec/2+1

#   Ensure that cb > vb
    if vb > cb:
        vb,cb = cb,vb
    elif vb == cb:
        cb = cb + 1
    if vb > nelec/2:
        print( "Warning: dealing with two empty bands: %3i  %3i" % (vb,cb))
    if cb < nelec/2:
        print( "Warning: dealing with two occupied bands: %3i  %3i" % (vb,cb))

    if vb == nelec/2 and cb == (nelec/2+1):
        print( "Calculate the minimal valence-conduction gap")
    if opts.kav:
        kavgap = vasp_anal_get_kavgap(band_struct,vb,cb,fix_k=opts.fix_k,inv=opts.flag_inv,debug=opts.debug)
    else:
        vasp_anal_get_gap(band_struct,vb,cb,debug=opts.debug)


# =====================================================

if __name__ == "__main__":
    Main(sys.argv)

