#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pa_init_input.py
# Creation Date : 18-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

from pa_classes import abinit_input_files
from argparse import ArgumentParser
import sys

def abinit_init_input(casename, formula, xc_type, pp_type, abinit_address, nproc, mpitype):
    
    a_input = abinit_input_files(casename, formula, xc_type, pp_type)
    a_input.set_abinit_run(abinit_address, nproc, mpitype)

    return a_input

# ====================================================

def Main(ArgList):

    description = '''Initialize the ABINIT input files and a running script'''
    
    parser = ArgumentParser(description=description)
    
    parser.add_argument(dest="casename", metavar="CASE", type=str, help="Casename for the calculation")
    parser.add_argument("-f", dest="formula", type=str, default='', help="Formula of the chemical system, e.g. Si1O2")
    parser.add_argument('--xc', dest='xc_type', type=str, default='PBE', help="type of XC functional. Default PBE")
    parser.add_argument('-n', dest='nproc', type=int, default=1, help="The number of processors used to run ABINIT. Default 1, i.e. serial calculation")
    parser.add_argument('--pp', dest='pp_type', type=str, default='paw', help="type of pseudopotential. Default PAW.")
    parser.add_argument('--mpi', dest='mpitype', type=str, default='mpiexec', help="Type of MPI executive. mpirun|mpiexec|yhrun")
    parser.add_argument('--ae', dest='abinit_address', default='abinit', type=str, help="The path of ABINIT executive, default by searching PATH.")
    
    # initialize options as 'opts'
    opts = parser.parse_args()

    xc_type = opts.xc_type.lower()
    pp_type = opts.pp_type.lower()
    
    a_input = abinit_init_input(opts.casename, opts.formula, xc_type, pp_type, \
                                opts.abinit_address, opts.nproc, opts.mpitype)

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

