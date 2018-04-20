#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pa_init_input.py
# Creation Date : 18-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from pa_classes import abinit_input_files
import sys



def abinit_init_input(casename, formula, xc_type, pp_type, **kwargs):
    
    a_input = abinit_input_files(opts.casename, opts.formula, xc_type, pp_type)

    return a_input

# ====================================================

def Main(ArgList):

    # add from argparse import ArgumentParser at head
    
    parser = ArgumentParser(description='')
    
    parser.add_argument("casename", metavar="CASE", type=str, help="Casename for the calculation")
    parser.add_argument("formula", metavar="FORMULA", type=str, help="Formula of the chemical system, e.g. Si1O2")
    parser.add_argument('--xc', dest='xc_type', type=str, default='PBE', help="type of XC functional. Default PBE")
    parser.add_argument('--pp', dest='pp_type', type=str, default='paw', help="type of pseudopotential. Default PAW.")
    
    # initialize options as 'opts'
    opts = parser.parse_args()

    xc_type = opts.xc_type.lower()
    pp_type = opts.pp_type.lower()
    
    a_input = abinit_init_input(opts.casename, opts.formula, xc_type, pp_type)


# ==============================

if __name__ == "__main__":
    Main(sys.argv)

