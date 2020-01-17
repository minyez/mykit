#!/usr/bin/env python3

import subprocess as sp
import os
from math import sqrt
from mykit.core.utils import run_cmd

def run_vasp(vasp_cmd=None, nproc=1, mpi=None, stdout=None, stderr=None):
    """perfrom vasp calculation
    """
    if nproc > 1 and mpi is None:
        mpi = "mpirun"
    if nproc < 0:
        mpi = None
    
    # search vasp or vasp_std if vasp_cmd is set to None
    if vasp_cmd is None:
        try:
            vasp_cmd = str(sp.check_output(["which", "vasp"], stderr=sp.STDOUT)).split('\n')[0]
        except sp.CalledProcessError:
            try:
                vasp_cmd = str(sp.check_output(["which", "vasp_std"], stderr=sp.STDOUT)).split('\n')[0]
            except sp.CalledProcessError:
                raise ValueError("Path for vasp or vasp_std not found. Need manual setting")

    if mpi is None:
        cmd = [vasp_cmd,]
    else:
        cmd = [mpi, "-np", str(nproc), vasp_cmd]
    
    run_cmd(cmd, fout=stdout, ferr=stderr)

def get_npar(nproc):
    """get NPAR parameter from number of processors
    """
    npar = 1
    nsqrt = int(sqrt(nproc)) + 1
    while nsqrt > 0:
        if nsqrt % 2 == 0 and nproc % nsqrt == 0:
            npar = nsqrt
            break
        nsqrt -= 1
    return npar
