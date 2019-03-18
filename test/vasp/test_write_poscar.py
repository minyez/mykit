#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from mykit.vasp.poscar import Poscar

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: %s poscar_file [...]" % sys.argv[0])
        sys.exit(-1)
    poscarFiles = sys.argv[1:]
    for pf in poscarFiles:
        pc = Poscar.read_from_file(pf)
        pc.write(pf+'_output')
