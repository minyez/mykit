#!/usr/bin/env python3
'''Prepare simple input files for a particular task.

If POSCAR exists, POTCAR and k-points will also be generated.
'''

import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from mykit.core.kmesh import kpath_encoder
from mykit.core.log import verbose
from mykit.core.symmetry import special_kpoints
from mykit.vasp.incar import incar
from mykit.vasp.kpoints import kpoints
from mykit.vasp.poscar import poscar
from mykit.vasp.potcar import potcar_search


def pv_simple_input():

    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", dest='poscar', default="POSCAR", \
        help="input file for KPOINTS and POTCAR generation. POSCAR as default")
    parser.add_argument("-e", dest='encut', type=int, default=None, \
        help="planewave cutoff.")
    parser.add_argument("-n", dest='nproc', type=int, default=1, \
        help="number of processors")
    parser.add_argument("-x", dest='xc', type=str, default=None, \
        help="type of XC functional for input orbitals, None for LEXCH in POTCAR")
    parser.add_argument("-l", dest='klen', type=int, default=0, \
        help="k-mesh density control, i.e. a*k. Negative for not generating KPOINTS")
    parser.add_argument("--spin", dest='ispin', type=int, default=1, \
        help="spin-polarization. 1 for nsp and 2 for sp.")
    parser.add_argument("--task", dest="task", \
        choices=["scf", "dos", "band", "gw", "opt"], \
        default="scf", \
        help="type of task for the input. Default scf.")
    parser.add_argument("-f", dest='overwrite', action="store_true", \
        help="flag to overwrite files, i.e. KPOINTS, POTCAR and INCAR")

    opts  = parser.parse_args()
    klen  = opts.klen
    klenDe = {"dos":40, "band": 15, "gw": 10}
    klenScf = 25

    if opts.overwrite:
        verbose.print_cm_warn("Overwrite switch on.")
    if os.path.isfile(opts.poscar):
        pos = poscar.read_from_file(opts.poscar)
        if klen >= 0 and (not os.path.isfile('KPOINTS') or opts.overwrite):
            if klen==0:
                klen = klenDe.get(opts.task, klenScf)
            if opts.task == "band":
                kpaths = special_kpoints.get_kpaths_from_cell(pos)
                if kpaths is None:
                    verbose.print_cm_warn("No predefined kpath available. Skip writing band KPOINTS.")
                else:
                    for i, kpath in enumerate(kpaths):
                        kpathStr = kpath_encoder(kpath["symbols"])
                        kp = kpoints(kmode="L", kdense=klen, kpath=kpath, \
                            comment="K-point path {}".format(kpathStr))
                        kp.write(pathKpoints="KPOINTS_band_{}".format(i))
            else:
                nks = [int(klen/x) for x in pos.alen]
                kp = kpoints(kmode="G", kgrid=nks)
                kp.write()
        if not os.path.isfile('POTCAR') or opts.overwrite:
            pts = potcar_search(*pos.atomTypes)
            pts.export(xc=opts.xc)

    # write INCAR
    if not os.path.isfile('INCAR') or opts.overwrite:
        ic = incar.minimal_incar(opts.task, xc=opts.xc, nproc=opts.nproc, \
            ISPIN=opts.ispin, ENCUT=opts.encut, comment="Quick {} input by mykit".format(opts.task))
        ic.write()
        


if __name__ == "__main__":
    pv_simple_input()
