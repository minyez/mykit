#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser

def pc_lattconv():
    """
    """
    parser = ArgumentParser(description=pc_lattconv.__doc__)
    group1 = parser.add_mutually_exclusive_group()
    
    parser.add_argument("pos", type=int, default=0, \
        help="template positional argument")
    parser.add_argument("-t", dest="test", type=int, default=0, \
        help="template optional argument")
    group1.add_argument("-x", dest="exclu", type=int, default=0, \
        help="template exclusive argument")
    parser.add_argument('-p', dest='paras', default=None, nargs='+', \
        help="optional argument with undetermined number of paras")
    
    # initialize options as 'opts'
    args = parser.parse_args()


if __name__ == "__main__":
    pc_lattconv()
