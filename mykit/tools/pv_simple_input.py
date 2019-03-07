#!/usr/bin/env python

from __future__ import print_function
from pv_calc_utils import vasp_io_get_NPAR, \
                          vasp_write_kpoints_basic, \
                          vasp_write_incar_minimal_elec
from pv_classes import vasp_read_poscar
from argparse import ArgumentParser
import os
import sys

# =====================================================

def Main(ArgList):

    description = '''
    Prepare simple INCAR file and k-points (if POSCAR exists)
    '''

# =================== Parser ==========================

    parser = ArgumentParser(description=description)
    parser.add_argument("-i",dest='poscar',default="POSCAR",help="Input file for KPOINTS generation. POSCAR as default")
    parser.add_argument("-e",dest='encut',type=int,default=0,help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors ")
    parser.add_argument("-x",dest='tag_xc',default=None,help="type of XC functional for input orbitals, \
                        None for LEXCH in POTCAR")
    parser.add_argument("-l",dest='klen',type=int,default=0,help="K-mesh density control, i.e. a*k. Negative for not generating KPOINTS")
    parser.add_argument("--spin",dest='ispin',type=int,default=1,help="Spin-polarization. 1 for nsp and 2 for sp.")

    opts  = parser.parse_args()
    poscar = opts.poscar
    klen  = opts.klen
    tag_xc= opts.tag_xc
    ispin = opts.ispin
    nproc = opts.nproc
    encut = opts.encut
    npar  = vasp_io_get_NPAR(nproc)

    print(" ============ pv_simple_input.py ============")

    if os.path.exists(poscar) and (klen>=0):
        print(" Writing KPOINTS...(check odd/even yourself!)")
        poscar = vasp_read_poscar(poscar)
        if klen==0:
            print("  - Warning: KLEN is not specified. Use default value: 30")
            print("  - Warning: you need to specify it yourself for metallic system")
            klen = 30
        nks = [int(klen/x) for x in poscar.lenlat]
        vasp_write_kpoints_basic(nks,'G')

    elif not os.path.exists(poscar):
        print(" - Error: need POSCAR file to generate KPOINTS. Exit")
        sys.exit(1)

    with open('INCAR','w') as incar:
        print(" Writing INCAR...")
        vasp_write_incar_minimal_elec(incar,tag_xc,encut=encut,npar=npar,spin=ispin)

# =====================================================

if __name__ == "__main__":
    Main(sys.argv)
