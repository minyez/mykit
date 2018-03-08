#!/usr/bin/env python

# ==============================
#     File Name : pv_addpot.py
# Creation Date : 16-10-2017
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : create POTCAR from potcar directory
# ==============================

import sys
import os
import subprocess as sp
from pv_classes import vasp_read_poscar
from argparse import ArgumentParser
import fnmatch

# ==============================

def Main(ArgList):

    #potdir = '/home/stevezhang/software/vasppot-5.2/'
    print " ============ pv_addpot.py ============"
    potdir = sp.check_output('echo $VASPPOT',shell=True).split()
    if len(potdir)==0:
        print "The environment variable $VASPPOT is not set. Exit."
        sys.exit(1)
    elif len(potdir)==1:
        potdir = potdir[0]

    if not os.path.exists(potdir):
        print "Invalid POTCAR directory"
        sys.exit(1)

    description = ' Create POTCAR from potcar directory: %s' % potdir
    parser = ArgumentParser(description=description)
    parser.add_argument('-e',dest='elements',default=None,nargs='+',help="Elements name, e.g. Cu, Fe_sv, H_GW, etc")
    parser.add_argument("-x",dest='xc',help="XC functional to generate PAW. LDA or PBE (default).",default='PBE')
    parser.add_argument("-c",dest='check',help="flag for check mode, checking available POTCARs",action='store_true')
    parser.add_argument('-i',dest='posin',default='POSCAR',help="POSCAR file for checking mode")
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')

    opts = parser.parse_args()
    if potdir.endswith('5.2') or potdir.endswith('5.3'):
    # nomenclature suits for vasp-5.2 and vasp-5.3
        potdir_xc = potdir+'/paw_'+opts.xc.upper()+'/'
    elif potdir.endswith('5.4'):
    # nomenclature suits for vasp-5.4
        potdir_xc = potdir+'/potpaw_'+opts.xc.upper()+'/'

    posin = opts.posin

    print " Finding POTCARs from: %s" % potdir_xc

    # check if the XC directory exists    
    if not os.path.exists(potdir_xc):
        print "POTCAR directory does not exist. Exit."
        sys.exit(1)
    if opts.elements is not None:
        ele = opts.elements
    else:
        poscar=vasp_read_poscar(posin)
        ele = poscar.atom_type

    if opts.debug:
        print potdir_xc
        print ele
        print opts.xc

    if opts.check:
        for iele in ele:
            print "%s avail:" % iele
            for idir in os.listdir(potdir_xc):
                if fnmatch.fnmatch(idir,iele):
                    print " - " + idir
                if fnmatch.fnmatch(idir,iele+"_*"):
                    print " - " + idir
        return

    with open('POTCAR','w') as fpot:
        for i in xrange(len(ele)):
            e_potcar = potdir_xc+ele[i]+'/POTCAR'
            if opts.debug:
                print e_potcar
            if os.path.exists(e_potcar):
                with open(e_potcar,'r') as fpotcar:
                    lines = fpotcar.readlines()
                    for line in lines:
                        fpot.write(line)
                print " - POTCAR of %8s found.     Add" % ele[i]
            else:
                print " - POTCAR of %8s not found. Skip" % ele[i]
    return

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

