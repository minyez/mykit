#!/usr/bin/env python

from pv_calc_utils import *
from pv_classes import vasp_read_poscar
import os

# =====================================================

def Main(ArgList):

    description = '''
    Prepare simple INCAR file and k-points (if POSCAR exists)
    '''

# =================== Parser ==========================

    parser = ArgumentParser(description=description)
    parser.add_argument("-i",dest='posin',default="POSCAR",help="Input file for KPOINTS generation. POSCAR as default")
    parser.add_argument("-e",dest='encut',type=int,default=0,help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors ")
    parser.add_argument("-x",dest='tag_xc',default=None,help="type of XC functional for input orbitals, \
                        None for LEXCH in POTCAR")
    parser.add_argument("-l",dest='klen',type=int,default=0,help="K-mesh density control, i.e. a*k. Negative for not generating KPOINTS")
    parser.add_argument("--spin",dest='ispin',type=int,default=1,help="Spin-polarization. 1 for nsp and 2 for sp.")

    opts  = parser.parse_args()
    posin = opts.posin
    klen  = opts.klen
    tag_xc= opts.tag_xc
    ispin = opts.ispin
    nproc = opts.nproc
    encut = opts.encut
    npar  = vasp_io_get_NPAR(nproc)

    print " ============ pv_simple_input.py ============"

    if os.path.exists(posin) and (klen>=0):
        print " Writing KPOINTS...(check odd/even yourself!)"
        poscar = vasp_read_poscar(posin)
        if klen==0:
            print "  - Warning: KLEN is not specified. Use default value: 30"
            print "  - Warning: you need to specify it yourself for metallic system"
            klen = 30
        nks = [int(klen/x) for x in poscar.lenlat]
        vasp_write_kpoints_basic(nks,'G')

    elif not os.path.exists(posin):
        print " - Error: need POSCAR file to generate KPOINTS. Exit"
        sys.exit(1)

    with open('INCAR','w') as incar:
        print " Writing INCAR..."
        vasp_write_incar_minimal_elec(incar,tag_xc,encut=encut,npar=npar,spin=ispin)

# =====================================================

if __name__ == "__main__":
    Main(sys.argv)
