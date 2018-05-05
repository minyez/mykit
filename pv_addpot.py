#!/usr/bin/env python

# ==============================
#     File Name : pv_addpot.py
# Creation Date : 16-10-2017
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : create POTCAR from potcar directory
# ==============================

from __future__ import print_function
import sys
import os
import subprocess as sp
from pv_classes import vasp_read_poscar
from argparse import ArgumentParser
import fnmatch


def vasp_getpotdir_xc(xc='PBE', verbose=False):

    potdir = sp.check_output('echo $VASPPOT',shell=True).split()

    if len(potdir)==0:
        if verbose:
            print("The environment variable $VASPPOT is not set. Exit.")
        sys.exit(1)
    elif len(potdir)==1:
        potdir = potdir[0]

    if not os.path.exists(potdir):
        if verbose:
            print("Invalid POTCAR directory")
        sys.exit(1)

    if potdir.endswith('5.2') or potdir.endswith('5.3'):
    # nomenclature suits for vasp-5.2 and vasp-5.3
        if xc == 'LDA':
            potdir_xc = potdir + '/paw_LDA/'
        else:
            potdir_xc = potdir + '/paw_PBE/'
    elif potdir.endswith('5.4'):
    # nomenclature suits for vasp-5.4
        if xc == 'LDA':
            potdir_xc = potdir + '/potpaw_LDA/'
        else: 
            potdir_xc = potdir + '/potpaw_PBE/'

    if not os.path.exists(potdir_xc):
        if verbose:
            print("POTCAR directory does not exist. Exit.")
        sys.exit(1)

    return potdir_xc


def vasp_addpot(xc, ele, backup=False, verbose=False):

    potdir_xc = vasp_getpotdir_xc(xc, verbose=False)

    if os.path.exists('POTCAR') and backup:
        os.rename('POTCAR','POTCAR.old')

    with open('POTCAR','w') as fpot:
        for i in xrange(len(ele)):
            e_potcar = potdir_xc+ele[i]+'/POTCAR'
            if os.path.exists(e_potcar):
                with open(e_potcar,'r') as fpotcar:
                    lines = fpotcar.readlines()
                    for line in lines:
                        fpot.write(line)
                if verbose:
                    print(" - POTCAR of %8s found.     Add" % ele[i])
            else:
                if verbose:
                    print(" - POTCAR of %8s not found. Skip" % ele[i])

# ==============================

def Main(ArgList):

    #potdir = '/home/stevezhang/software/vasppot-5.2/'

    description = ' Create POTCAR from potcar directory defined as environment variable VASPPOT, e.g. /home/stevezhang/software/vasppot-5.2/ with paw_LDA and paw_PBE as subdirectories.'

    parser = ArgumentParser(description=description)
    parser.add_argument('-e',dest='elements',default=None,nargs='+',help="Elements name, e.g. Cu, Fe_sv, H_GW, etc")
    parser.add_argument("-x",dest='xc',help="XC functional to generate PAW. LDA or PBE (default).",default='PBE')
    parser.add_argument("-c",dest='check',help="flag for check mode, checking available POTCARs",action='store_true')
    parser.add_argument('-i',dest='posin',default='POSCAR',help="POSCAR file for checking mode")
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')

    opts = parser.parse_args()

    print(" ============ pv_addpot.py ============")
    #potdir = sp.check_output('echo $VASPPOT',shell=True).split()
    #if len(potdir)==0:
    #    print("The environment variable $VASPPOT is not set. Exit.")
    #    sys.exit(1)
    #elif len(potdir)==1:
    #    potdir = potdir[0]

    #if not os.path.exists(potdir):
    #    print("Invalid POTCAR directory")
    #    sys.exit(1)

    #if potdir.endswith('5.2') or potdir.endswith('5.3'):
    ## nomenclature suits for vasp-5.2 and vasp-5.3
    #    potdir_xc = potdir+'/paw_'+opts.xc.upper()+'/'
    #elif potdir.endswith('5.4'):
    ## nomenclature suits for vasp-5.4
    #    potdir_xc = potdir+'/potpaw_'+opts.xc.upper()+'/'

    # check if the XC directory exists    
    if opts.elements is not None:
        ele = opts.elements
    else:
        posin = opts.posin
        poscar=vasp_read_poscar(posin)
        ele = poscar.atom_type

    if opts.debug:
        print(opts.xc)
        print(ele)

    if opts.check:
        potdir_xc = vasp_getpotdir_xc(opts.xc, verbose=True)
        print(" Finding POTCARs from: %s" % potdir_xc)
        for iele in ele:
            print("%s avail:" % iele)
            for idir in os.listdir(potdir_xc):
                if fnmatch.fnmatch(idir,iele):
                    print(" - " + idir)
                if fnmatch.fnmatch(idir,iele+"_*"):
                    print(" - " + idir)
    else:
        vasp_addpot(opts.xc, ele, backup=True, verbose=True)

    return

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

