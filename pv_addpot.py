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
import fnmatch
import subprocess as sp
from argparse import ArgumentParser
from pv_classes import vasp_read_poscar


def vasp_getpotdir_xc(xc='PBE', verbose=False):

    potdir = sp.check_output('echo $VASPPOT', shell=True).split()

    if len(potdir) == 0:
        if verbose:
            print("The environment variable $VASPPOT is not set. Exit.")
        sys.exit(1)
    
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
            potdir_xc = potdir + '/potpaw'
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
        os.rename('POTCAR', 'POTCAR.old')

    with open('POTCAR', 'w') as fpot:
        for i in range(len(ele)):
            e_potcar = potdir_xc + ele[i] + '/POTCAR'
            if os.path.exists(e_potcar):
                with open(e_potcar, 'r') as fpotcar:
                    lines = fpotcar.readlines()
                    for line in lines:
                        fpot.write(line)
                if verbose:
                    print(" - POTCAR of %8s found.     Add" % ele[i])
            else:
                if verbose:
                    print(" - POTCAR of %8s not found. Skip" % ele[i])

def __get_enmax(potcarpath):
    '''
    '''
    return float(sp.check_output(['grep', 'ENMAX', potcarpath]).split()[2][:-1])

def __get_enmin(potcarpath):
    '''
    '''
    return float(sp.check_output(['grep', 'ENMIN', potcarpath]).split()[5])

# ==============================

def Main(ArgList):

    #potdir = '/home/stevezhang/software/vasppot-5.2/'

    '''
    Create POTCAR from potcar directory defined as environment variable VASPPOT,
    e.g. /home/stevezhang/software/vasppot-5.2/ with paw_LDA and paw_PBE as subdirectories.
    '''

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('-e', dest='elements', default=None, nargs='+', \
            help="Elements name, e.g. Cu, Fe_sv, H_GW, etc")
    parser.add_argument("-x", dest='xc', default='PBE', \
            help="XC functional to generate PAW. LDA or PBE (default).")
    parser.add_argument("-c", dest='check', action='store_true', \
            help="flag for check mode, checking available POTCARs")
    parser.add_argument('-i', dest='posin', default='POSCAR', \
            help="POSCAR file for checking mode")
    parser.add_argument("-D", dest='debug', action='store_true', \
            help="flag for debug mode")

    opts = parser.parse_args()

    print(" ============ pv_addpot.py ============")

    if opts.elements is not None:
        ele = opts.elements
    else:
        posin = opts.posin
        poscar = vasp_read_poscar(posin)
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
                if fnmatch.fnmatch(idir, iele) or fnmatch.fnmatch(idir, iele+"_*"):
                    enmax = __get_enmax(os.path.join(potdir_xc, idir + '/POTCAR'))
                    enmin = __get_enmin(os.path.join(potdir_xc, idir + '/POTCAR'))
                    print(" - %-12s (%.2f - %.2f eV)" % (idir, enmin, enmax))
    else:
        vasp_addpot(opts.xc, ele, backup=True, verbose=True)

    return

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

