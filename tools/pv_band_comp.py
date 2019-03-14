#!/usr/bin/env python
#
# a simple script to get atom-l component of band \psi_{nk} at each k-point

from pv_classes import vasp_read_xml
from argparse import ArgumentParser
import numpy as np
import os,sys

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Calculate the l-weighted gap, either minimal or k-averaged between two bands
    Also these two values between any occupied and empty band is possible
    by adding option: -v vb -c cb
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')
    parser.add_argument("-n",dest='ib',help="the index of the target band",type=int,default=0)
    parser.add_argument("-l",dest='ldecomp',help="show the l-decomposed band component",default=None)
    parser.add_argument("-s",dest='sdecomp',help="show the atomtype-decomposed band component",default=None)

    opts = parser.parse_args()
    dxml = vasp_read_xml()
    nelec   = dxml.nelec
    nkp     = dxml.nkp
    bandmax = dxml.nbands
    assert -1 <= ib <= bandmax
    if ib == 0:
        ib = nelec/2   #VBM
    elif ib == -1:
        ib = nelec/2 + 1   #CBM
    else:
        ib = opts.ib

# =====================================================


#   read transition
#    if len(vc) != 4:
#        print "Invalid transition."
#        sys.exit(1)
#    else:
#        try:
#            atomtype_v = int(vc[0])
#            orb_v = vc[1]
#            atomtype_c = int(vc[2])
#            orb_c = vc[3]
#        except:
#            print "Invalid transition."
#            sys.exit(1)

    vb_pw_index = dxml.pwave_index(orb_v)
    vb_at_index = dxml.atoms_index(atomtype_v)
    cb_pw_index = dxml.pwave_index(orb_c)
    cb_at_index = dxml.atoms_index(atomtype_c)

    if opts.debug:
        print atomtype_v,orb_v,vb_pw_index,vb_at_index
        print atomtype_c,orb_c,cb_pw_index,cb_at_index

    try:
        assert vb >= 0 and cb >= 0
    except AssertionError:
        print "Invalid band index"
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

    kgap = []
    for kp in xrange(dxml.nkp):
        vb_atom_l_weigh = dxml.sum_atom_l_comp(0,kp,vb-1,vb_at_index,vb_pw_index)
        cb_atom_l_weigh = dxml.sum_atom_l_comp(0,kp,cb-1,cb_at_index,cb_pw_index)
        gap = dxml.eigen[0][kp][cb-1] - dxml.eigen[0][kp][vb-1]
        if opts.flag_inv:
            gap = 1/gap
        if opts.debug:
            print gap,vb_atom_l_weigh,cb_atom_l_weigh
        kgap.append(gap*vb_atom_l_weigh*cb_atom_l_weigh)

    kavgap = np.inner(np.array(kgap),dxml.kpweigh)/np.sum(dxml.kpweigh)
    if opts.flag_inv:
        kavgap = 1.0E0/kavgap
#    print kavgap
#    print np.array(kgap)
#    print dxml.kpweigh
#    print sum(dxml.kpweigh)
#    print kavgap


#   Ensure that cb > vb
    if vb > cb:
        vb,cb = cb,vb
    elif vb == cb:
        cb = cb + 1
    if vb > nelec/2:
        print "Warning: dealing with two empty bands: %3i  %3i" % (vb,cb)
    if cb < nelec/2:
        print "Warning: dealing with two occupied bands: %3i  %3i" % (vb,cb)

    if vb == nelec/2 and cb == (nelec/2+1):
        print "Calculate the minimal valence-conduction gap."
    if opts.flag_inv:
        print "Calculate inverse of band gap."
    if opts.kav:
        print "k-averged gap between band %i and %i: "%(vb,cb) ,kavgap
        pass
    else:
#        vasp_anal_get_gap(band_struct,vb,cb,debug=opts.debug)
        pass


# =====================================================

if __name__ == "__main__":
    Main(sys.argv)

