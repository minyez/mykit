#!/usr/bin/env python3
# coding=utf-8
"""compute the wavefunction of Local Orbitals (LO) by reading energy.

The following files are required

- casename.struct
- casename.in1 or casename.in1c
- casename.vsp

spin-polarized case is not supported yet.
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
from mykit.core.log import Verbose
from mykit.wien2k.struct import Struct
from mykit.wien2k.inputs import In1
from mykit.wien2k.outputs import Vsp
from mykit.wien2k.utils import find_complex_file, outwin
from mykit.core.utils import conv_string

def pw_plot_apw():
    """main stream of pw_plot_apw
    """

    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-c', dest='casename', type=str, default=None, \
        help="casename. Use filename of first struct file in PWD as default")
    parser.add_argument('--in1', dest='in1', type=str, default=None, \
        help="in1 file. Overwrite that specified by casename")
    parser.add_argument('--ali', dest='ali', type=str, nargs="+", \
        help="Index of LO, in format of 'a:l:i', a for atom index, l for angular moment and i for index of exception in l channel")
    #parser.add_argument('--ur', dest='plotur', action='store_true', \
    #    help="Plot wavefunction raidial product, i.e. u*r instead of u")
    #parser.add_argument('--single', dest='single', action='store_true', \
    #    help="Plot only the wavefunction of energy E")

    opts = parser.parse_args()
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    pathStruct = None
    pathIn1 = None
    pathVsp = None
    cn = opts.casename
    if cn is not None:
        pathStruct = cn + '.struct'
        pathVsp = cn + '.vsp'
        pathIn1 = find_complex_file(cn, 'in1')
    if opts.in1 is not None:
        pathIn1 = opts.in1
    
    struct = Struct.read_from_file(pathStruct=pathStruct)
    in1 = In1.read_from_file(filePath=pathIn1)
    vsp = Vsp(pathVsp)
    ali = opts.ali
    if ali is None:
        ali = []

    for x in ali:
        ia, l, idl = conv_string(x, int, sep=":")
        els = in1.get_exceptions(ia)
        if els is None or l not in els:
            Verbose.print_cm_warn("atom %d l %d not available" % (ia, l))
            continue
        nel = len(els[l])
        if idl >= nel:
            Verbose.print_cm_warn("%d-th orbit in %d l-channel not available" % (idl, l))
            continue
        E = els[l][idl][0]
        if E == 0.3:
            E = in1.efermi - 0.2
        E /= 2.0
        r = struct.get_radial(ia)
        vr = vsp.get_vsp(ia)
        Z = struct.get_z(ia)
        urlarge, ursmall, u, du, nodes = outwin(r, vr, E, l, Z, normalize=True)
        print("iat = ", ia, "l = ", l, "El(Ha) = ", E, "nodes = ", nodes)
        print(*urlarge[0:5])
        print(*ursmall[0:5])
        # if not opts.plotur:
        urlarge = np.divide(urlarge, r)
        ax.plot(r, urlarge, marker='.', ls='-', lw=2)

    plt.show()

if __name__ == "__main__":
    pw_plot_apw()