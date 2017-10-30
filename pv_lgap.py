#!/usr/bin/env python
#
# a simple script to get fundamental gap and k-point-averaged band gap

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
    parser.add_argument("-k",dest='kav',help="flag for print k-averaged gap",action='store_true')
    parser.add_argument("-v",dest='vb',help="the index of valence band",type=int,default=0)
    parser.add_argument("-c",dest='cb',help="the index of conduction band",type=int,default=0)
#   fix_k not supported yet
#    parser.add_argument("-f",dest='fix_k',help="flag for k-point fixing mode for k-averaged gap. 0 for fixing VBM, 1 for fixing CBM ",type=int,default=-1)
    parser.add_argument("-i",dest='flag_inv',help="flag for the inverse of gap. Useful for analysis of polarization",action='store_true')
    parser.add_argument("-l",dest='vc',help="the atom-l-components over which the valence and the conduction band will be weighed, respectively. E.g. '1s2p' for the transition from type1atom-s orbital to type2atom-p orbital. m-decomposition is not supported yet.",default='0t0t')
    parser.add_argument("-s",dest='showcomp',help="flag for show the band component at each k-point",action='store_true')


    opts = parser.parse_args()
    vb = opts.vb
    cb = opts.cb
    vc = opts.vc.strip()
    dxml = vasp_read_xml()

# =====================================================

#    band_struct = vasp_anal_read_eigen(spinpolarzied=False,debug=opts.debug)
    nelec   = dxml.nelec
    nkp     = dxml.nkp
    bandmax = dxml.nbands

#   read transition
    if len(vc) != 4:
        print "Invalid transition."
        sys.exit(1)
    else:
        try:
            atomtype_v = int(vc[0])
            orb_v = vc[1]
            atomtype_c = int(vc[2])
            orb_c = vc[3]
        except:
            print "Invalid transition."
            sys.exit(1)

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
    if opts.showcomp:
        if atomtype_v == 0:
            typevb = 'All'
        if atomtype_c == 0:
            typecb = 'All'
        if atomtype_v != 0 and atomtype_c != 0:
            typevb = dxml.type[atomtype_v-1]
            typecb = dxml.type[atomtype_c-1]
        print "k-point         k-vec             Eg_direct      VB %i:%s-%s     CB %i:%s-%s" % (vb,typevb,orb_v,cb,typecb,orb_c)
    for kp in xrange(dxml.nkp):
    # currently only the first spin component
        vb_atom_l_weigh = dxml.sum_atom_l_comp(0,kp,vb-1,vb_at_index,vb_pw_index)
        cb_atom_l_weigh = dxml.sum_atom_l_comp(0,kp,cb-1,cb_at_index,cb_pw_index)
        gap = dxml.eigen[0][kp][cb-1] - dxml.eigen[0][kp][vb-1]
        if opts.debug:
            print gap,vb_atom_l_weigh,cb_atom_l_weigh
        if opts.showcomp:
            print "%7s  (%5.3f,%5.3f,%5.3f)      %9.4f        %9.4f       %9.4f" % \
                 (kp+1,dxml.kplist[kp][0],dxml.kplist[kp][1],dxml.kplist[kp][2],gap,vb_atom_l_weigh,cb_atom_l_weigh)
        if opts.flag_inv:
            gap = 1.0/gap
        kgap.append(gap*vb_atom_l_weigh*cb_atom_l_weigh)

    kavgap = np.inner(np.array(kgap),dxml.kpweigh)/np.sum(dxml.kpweigh)
    if opts.flag_inv:
        kavgap = 1.0E0/kavgap


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

