#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pv_addvac.py
# Creation Date : 30-10-2017
# Last Modified : Mon 30 Oct 2017 10:01:30 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose :
#
# ====================================================

import sys
from argparse import ArgumentParser
from pv_classes import vasp_read_poscar

def Main(ArgList):
    description = 'Change the thickness of vacuum of the slab model'
    parser = ArgumentParser(description=description)
    parser.add_argument("-d",dest='zdirt',type=int,default=3,help="the direction of non-periodic dimension")
    parser.add_argument("-v",dest='vac',type=float,default=0.0,\
                        help="the thichness of vacuum to add (negative to reduce). In Angstrom")
    parser.add_argument("-i",dest='input',default='POSCAR',help="input file. Default POSCAR")
    parser.add_argument("-o",dest='output',default='POSCAR_new',help="output file. Default POSCAR_new")

    opts = parser.parse_args()
    print " ============ pv_addvac.py ============"
    poscar = vasp_read_poscar(opts.input)
    iz = opts.zdirt - 1
    print "   Original z: %17.8f" % poscar.lattice[iz][iz]
    poscar.action_add_vacuum(opts.vac,opts.zdirt)
    print "   Changed z: %18.8f" % poscar.lattice[iz][iz]
    poscar.write_poscar(opts.output)


# ==============================

if __name__ == "__main__":
    Main(sys.argv)

