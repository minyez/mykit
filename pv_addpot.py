#!/usr/bin/env python

# ==============================

#     File Name : pv_addpot.py
# Creation Date : 16-10-2017
# Last Modified : Mon 16 Oct 2017 06:33:49 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : create POTCAR from potcar directory
#
# ==============================

import sys
import os
from argparse import ArgumentParser

# ==============================

def Main(ArgList):

    potdir = '/home/stevezhang/software/vasppot-5.2/'
    description = '''Create POTCAR from potcar directory: %s''' % potdir

    parser = ArgumentParser(description=description)
    parser.add_argument('elements',nargs='+',help="Elements name, e.g. Cu, Fe_sv, H_GW, etc")
    parser.add_argument("-x",dest='xc',help="XC functional to generate PAW. LDA or PBE (default).",default='PBE')
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')

    opts = parser.parse_args()
    potdir_xc = potdir+'paw_'+opts.xc.upper()+'/'
    ele = opts.elements
    if opts.debug:
        print potdir_xc
        print opts.elements
        print opts.xc

    if not os.path.exists(potdir_xc):
        print "POTCAR directory does not exist. Exit."
        sys.exit(1)
    else:
        print description
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
                    print "POTCAR of %8s found.     Add" % ele[i]
                else:
                    print "POTCAR of %8s not found. Skip" % ele[i]

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

