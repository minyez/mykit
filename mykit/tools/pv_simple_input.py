#!/usr/bin/env python3
'''Prepare simple input files for a particular task.

If POSCAR exists, POTCAR and k-points will also be generated.
'''

import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from mykit.vasp.incar import incar
from mykit.vasp.kpoints import kpoints
from mykit.vasp.poscar import poscar
from mykit.vasp.potcar import potcar_search


def pv_simple_input():


    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", dest='poscar', default="POSCAR", \
        help="Input file for KPOINTS generation. POSCAR as default")
    parser.add_argument("-e", dest='encut', type=int, default=None, \
        help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n", dest='nproc', type=int, default=1, \
        help="Number of processors ")
    parser.add_argument("-x", dest='xc', type=str, default=None, \
        help="type of XC functional for input orbitals, None for LEXCH in POTCAR")
    parser.add_argument("-l", dest='klen', type=int, default=0, \
        help="K-mesh density control, i.e. a*k. Negative for not generating KPOINTS")
    parser.add_argument("--spin", dest='ispin', type=int, default=1, \
        help="Spin-polarization. 1 for nsp and 2 for sp.")
    parser.add_argument("--task", dest="task", type=str, default="scf", \
        help="type of task for the input. Default scf (TODO)")

    opts  = parser.parse_args()
    klen  = opts.klen

    # try:
    #     assert os.path.isfile(opts.poscar)
    # except AssertionError:
    #     raise FileNotFoundError("need POSCAR file to generate KPOINTS")

    if os.path.exists(opts.poscar) and (klen>=0):
        # print(" Writing KPOINTS...(check odd/even yourself!)")
        pos = poscar.read_from_file(opts.poscar)
        if klen==0:
            # print("  - Warning: KLEN is not specified. Use default value: 30")
            # print("  - Warning: you need to specify it yourself for metallic system")
            klen = 30
        nks = [int(klen/x) for x in pos.alen]
        kp = kpoints(kmode="G", kgrid=nks)
        kp.write()
        pts = potcar_search(*pos.atomTypes)
        pts.export(xc=opts.xc)

    # write INCAR
    ic = incar.minimal_incar(opts.task, xc=opts.xc, nproc=opts.nproc, \
        ISPIN=opts.ispin, ENCUT=opts.encut, comment="Simple {} input by mykit".format(opts.task))
    ic.write()


if __name__ == "__main__":
    pv_simple_input()
