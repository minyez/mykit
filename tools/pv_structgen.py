#!/usr/bin/env python3
'''
Generate a number of structures with different volumes
for plotting energy-volume curve in VASP
'''

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
from mykit.vasp.poscar import Poscar
from mykit.core.utils import get_arith_prog

def pv_structgen():
    '''main stream
    '''

    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", dest='ipos', type=str, default='POSCAR', \
        help="input struct file, default POSCAR")
    parser.add_argument("-s", dest='stvra', type=float, default=0.92, \
        help="start volume ratio, default 0.92")
    parser.add_argument("-e", dest='edvra', type=float, default=1.08, \
        help="end volume ratio, default 1.08")
    parser.add_argument("-d", dest='interval', type=float, default=0.02, \
        help="interval of volume ratio, default 0.02")
    parser.add_argument("--dir", action="store_true", \
        help="save each POSCAR in directory with name 'V_x.xx'")

    opts = parser.parse_args()
    
    vras = get_arith_prog(opts.stvra, opts.edvra, interval=opts.interval)
    aras = np.power(vras, 1.0/3)

    pc = Poscar.read_from_file(opts.ipos) 

    if opts.dir:
        for i, x in enumerate(vras):
            dname = 'V_%4.2f' % x
            if os.path.isdir(dname):
                print("%s found. Remove old..." % dname)
                os.removedirs(dname)
            os.makedirs(dname)
            pc.write(os.path.join(dname, 'POSCAR'), scale=aras[i])
    else:
        for i, x in enumerate(vras):
            fname = 'POSCAR_%3.2f' % x
            pc.write(fname, scale=aras[i])
    

if __name__ == "__main__":
    pv_structgen()
