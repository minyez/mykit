#!/usr/bin/env python
# use setrmt_lapw and clminter to restart the wien2k calculation from a different RMT setting

from pw_anal_utils import Get_Casename
import sys,os
import subprocess as sp
from argparse import ArgumentParser
from shutil import copy2


def Main(ArgList):
    example = " Example: pw_resetrmt Cu:2.1 Cl:2.1."
    warning = ' Please input in the format of ``Element:RMT'' one by one.\n %s' %example
    description = ' Reset the RMT of case.struct.\n%s' % warning
    parser = ArgumentParser(description=description)

    parser.add_argument('elements',nargs='+',help="Elements name and its RMT, in the form of Ele:Rmt, separated by space")
    parser.add_argument("-f",dest="casename",help="Case name",default=None)
    parser.add_argument("-i",action='store_true',help="Flag for interpolate the charge density from results of old RMT")

    opts = parser.parse_args()

    if opts.casename is None:
        casename = Get_Casename()
        case_struct = casename+'.struct'
    else:
        casename = opts.casename
        case_struct = casename.strip()+'.struct'

    new_rmt = opts.elements

    if len(new_rmt) == 0:
        print "Current case: " + casename
        print warning
        sys.exit(1)
    if not os.path.isfile(case_struct):
        print " %s.struct file not found." % casename
        sys.exit(1)

    for i in xrange(len(new_rmt)):
        new_rmt[i] = new_rmt[i].strip()
        while new_rmt[i].startswith(','):
            new_rmt[i] = new_rmt[i][1:].strip()
        while new_rmt[i].endswith(','):
            new_rmt[i] = new_rmt[i][:-1].strip()

    new_rmt_str = ','.join(new_rmt)
    sp.call("setrmt_lapw "+casename+" -a "+new_rmt_str,shell=True)
    copy2(case_struct+'_setrmt',case_struct)
    os.remove(case_struct+'_setrmt')

    # interpolate charge density
    # for spin-polarized case, also need x clminter -up/-dn
    if opts.i:
        copy2(case_struct+'_setrmt',case_struct+'_new')
        sp.call("x clminter",shell=True)
        copy2(case_struct+'_new',case_struct)
        copy2(casename+'.clmsum_new',casename+'.clmsum')

# ====================================================

if __name__ == "__main__":
    Main(sys.argv)
