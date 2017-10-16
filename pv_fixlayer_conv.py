#!/usr/bin/env python

from pv_calc_utils import *
from pv_anal_utils import vasp_read_poscar
from shutil import copy2
from argparse import ArgumentParser
import subprocess as sp
import os

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Test the convergence of slab calculation w.r.t fixed layers
    Fix the bottom side. Asymmetric slab only. Need further implementation'''

    parser = ArgumentParser(description=description)
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors ")
#    parser.add_argument("-r",dest="range",type=int,nargs=2,help="The range of fixed layers")
    parser.add_argument("-s",dest='sym',help="Flag for symmetric fixing",action="store_true")
    parser.add_argument("-v",dest='vasp_path',default="vasp",help="Path of vasp executive")
    parser.add_argument("-D",dest='debug',help="Debug mode",action="store_true")

    opts = parser.parse_args()
#    npar = vasp_io_get_NPAR(opts.nproc)

    case_pos = vasp_read_poscar()
    natoms = case_pos.natoms
    fix_list = []
    rev_list = reversed(fix_list)

    if opts.debug:
        print "Total number of atoms: %d" % natoms

    for i in rev_list:
        os.chdir(str(i)+'fixed')
        if not i == 10:
            copy2('../'+str(i+1)+'fixed/WAVECAR','WAVECAR')
            copy2('../'+str(i+1)+'fixed/CONTCAR','CONTCAR')
            opt_num = int(sp.check_output("grep -c 'T  T  T' POSCAR_new",shell=True))
            fix_num = natoms - opt_num
            if opts.debug:
                print fix_num,type(fix_num)
            if not opts.sym:
                sp.call("pv_fix_slab.py -f CONTCAR -o POSCAR -n %d" % fix_num,shell=True)
            else:
                sp.call("pv_fix_slab.py -f CONTCAR -o POSCAR -n %d -s" % int(fix_num)/2,shell=True)
        vasp_path, vasp_cmd = vasp_vaspcmd_zmy(opts.nproc,vasp_path=opts.vasp_path)
        if opts.debug:
            print opt_num, fix_num, vasp_cmd
#        vasp_vasprun_zmy(vasp_cmd,'out','error')
        os.chdir('..')

# ====================================================

if __name__ == "__main__":
    Main(sys.argv)
