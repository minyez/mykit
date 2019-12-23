#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Workflow of WIEN2k optic calculation
"""

from argparse import ArgumentParser


def pw_optic():
    """Perform WIEN2k optic calculation
    """
    parser = ArgumentParser(description=pw_optic.__doc__)
    
    parser.add_argument("-t", dest="test", type=int, default=0, \
        help="template optional argument")
    parser.add_argument('-p', dest='paras', default=None, nargs='+', \
        help="optional argument with undetermined number of paras")
    
    # initialize options as 'opts'
    args = parser.parse_args()



if __name__ == "__main__":
    pw_optic()
