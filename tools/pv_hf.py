#!/usr/bin/env python3
# coding=utf-8

from __future__ import print_function
from mykit.vasp.incar import Incar
from mykit.vasp.utils import get_npar
from argparse import ArgumentParser
import os, sys, re
from shutil import copy2

def pv_write_hf_sequence(vasp_cmd, nproc, mpi, seq_out='hf_calc.py'):
    '''Write the python control file for a fast-SCF-converging HF calculation.
    '''

    # check if the three INCARs exist
    incar_list = ['INCAR_pbe', 'INCAR_nkred', 'INCAR_hf']
    for i in incar_list:
        if not os.path.isfile(i):
            raise ValueError("Error: specified INCAR not found, %s " % i)

    _template = """#!/usr/bin/env python3
import os
import subprocess as sp
from shutil import copy2
from mykit.vasp.utils import run_vasp
from mykit.core.utils import io_cleandir

vasp_cmd = '{vasp_cmd}'
nproc = {nproc}
mpi = '{mpi}'

ifiles = ["POSCAR", "POTCAR", "KPOINTS"]
for i in ifiles:
    if not os.path.exists(i):
        raise IOError("%s not found" % i)
# copy input files except INCAR
workdirs = ["1_PBE_preconv", "2_HF_coarse", "3_HF_normal"]
for w in workdirs:
    io_cleandir(w)
    for i in ifiles:
        copy2(i, w)
copy2('INCAR_pbe', workdirs[0]+"/INCAR")
copy2('INCAR_nkred', workdirs[1]+"/INCAR")
copy2('INCAR_hf', workdirs[2]+"/INCAR")

# PBE preconverge
os.chdir(workdirs[0])
run_vasp(vasp_cmd, nproc, mpi, "out", "error")

# coarse HF
os.chdir('../'+workdirs[1])
try:
    copy2("../"+workdirs[0]+"/WAVECAR", "WAVECAR")
except IOError:
    raise IOError("WARNING: WAVECAR in the PBE step not found")
run_vasp(vasp_cmd, nproc, mpi, "out", "error")

# standard HF
os.chdir('../'+workdirs[2])
try:
    copy2("../"+workdirs[1]+"/WAVECAR", "WAVECAR")
except IOError:
    raise IOError("WARNING: WAVECAR in the coarse HF step not found")
run_vasp(vasp_cmd, nproc, mpi, "out", "error")

# finish
os.chdir('../')
"""
    # write the python control file
    with open(seq_out, 'w') as h_out:
        h_out.write(_template.format(vasp_cmd=vasp_cmd, nproc=nproc, mpi=mpi))


def pv_prepare_hf_calculation():

    description = '''Prepare the input files and python executive script for a hybrid functional caluclation (PBE0, HSE06, HF).'''

    parser = ArgumentParser(description=description)
    group1 = parser.add_mutually_exclusive_group()

    parser.add_argument("-e", dest='encut', type=int, default=0, \
        help="planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n", dest='nproc', type=int, default=1, \
        help="number of processors")
    parser.add_argument("--mpi", dest='mpi', type=str, default="mpirun", \
        help="mpirun type")
    parser.add_argument("-x", dest='tag_xc', default='HSE06', \
        help="yype of hybrid functional. HSE06|PBE0|HF")
    parser.add_argument("-v", dest='vasp_cmd', default="vasp", \
        help="vasp executive")
    group1.add_argument("--nkred", dest='nkred', type=int, default=1, \
        help="NKRED set in the coarse calculation of HF ")
    group1.add_argument("--nkredxyz",dest='nkredxyz', \
        help="NKREDX/Y/Z set in the coarse calculation of HF, [X,Y,Z]")
    parser.add_argument("--spin", dest='ispin', type=int, default=1, \
        help="Spin-polarization. 1 for nsp and 2 for sp.")

    opts  = parser.parse_args()
    ispin = opts.ispin
    nproc = opts.nproc
    encut = opts.encut

    npar  = get_npar(nproc)

    # generate primitive INCAR files
    ic1 = Incar.minimal_incar("scf",
                              xc="PBE",
                              nproc=nproc,
                              ISPIN=ispin,
                              ENCUT=encut,
                              npar=npar,
                              comment="PBE preconv for hybrid")
    ic2 = Incar.minimal_incar("scf",
                              xc=opts.tag_xc,
                              nproc=nproc,
                              ISPIN=ispin,
                              ENCUT=encut,
                              npar=1,
                              PRECFOCK='FAST',
                              ISTART=1,
                              comment="hybrid calculation with nkred")
    ic3 = Incar.minimal_incar("scf",
                              xc=opts.tag_xc,
                              nproc=nproc,
                              ISPIN=ispin,
                              ENCUT=encut,
                              npar=1,
                              ISTART=1,
                              comment="hybrid calculation")
    if ic2["ALGO"] == "ALL":
        ic2["ISMEAR"] = 0
        ic2["SIGMA"] = 0.05
    if ic3["ALGO"] == "ALL":
        ic3["ISMEAR"] = 0
        ic3["SIGMA"] = 0.05


    # set coarse HF calculation with NKRED or NKREDX/Y/Z
    if opts.nkred != 1:
        ic2['NKRED'] = opts.nkred
    if opts.nkredxyz:
        try:
            nkredxyz = [x.strip() for x in re.split(r'[,\[\]]',opts.nkredxyz)]
            nkredxyz = [int(x) for x in ' '.join(nkredxyz).split()]

            nkred_opt = ['NKREDX','NKREDY','NKREDZ']
            for i in range(3):
                if nkredxyz[i] !=1 :
                    ic2[nkred_opt[i]] = nkredxyz[i]
        except:
            print(" WARNING: Unsupported input of NKREDX/Y/Z. Pass")

    ic1.write('INCAR_pbe')
    ic2.write('INCAR_nkred')
    ic3.write('INCAR_hf')

    pv_write_hf_sequence(nproc=opts.nproc, vasp_cmd=opts.vasp_cmd, mpi=opts.mpi)

# =====================================================

if __name__ == "__main__":
    pv_prepare_hf_calculation()

