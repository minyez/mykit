#!/usr/bin/env python3
'''Merge parallel files of wien2k by joinvec
'''

from mykit.wien2k.utils import get_casename
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import subprocess as sp

def pw_merge_ene():

    parser = ArgumentParser(description=__doc__,
        formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--so', dest='soc', action="store_true", \
        help="trigger SO mode")
    opts = parser.parse_args()

    cmds = ['x', 'joinvec']
    if opts.soc:
        cmds.append('-so')
    print("Running command:", *cmds)
    sp.check_output(cmds)


if __name__ == "__main__":
    pw_merge_ene()

