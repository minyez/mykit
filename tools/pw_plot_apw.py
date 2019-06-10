#!/usr/bin/env python3
# coding=utf-8
"""compute the wavefunction of Local Orbitals (LO) by reading energy.

The following files are required

- casename.struct
- casename.in1 or casename.in1c
- casename.vsp

spin-polarized case is not supported yet.
"""
from copy import deepcopy
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
from mykit.core.log import Verbose
from mykit.wien2k.struct import Struct
from mykit.wien2k.inputs import In1
from mykit.wien2k.outputs import Vsp
from mykit.wien2k.utils import find_complex_file, solve_rde, logr_int, solve_ene_deriv
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
        help="Index of LOs to plot, in format of 'a:l:i', a for atom index, l for angular moment and i for index of exception in l channel")
    parser.add_argument('--uonly', dest='uonly', action='store_true', \
        help="plot u of energy E only, instead of APW, lo or LO")
    parser.add_argument('-D', dest='debug', action='store_true', \
        help="debug mode")

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
    # remove duplicates
    ali = list(set(ali))
    requestial = []
    
    vrs = {}
    rs = {}
    Zs = {}
    los = {}
    apw = {}

    # retrive the linearization energy of LOs
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
        requestial.append((ia, l))
        E = els[l][idl][0]
        if E == 0.3:
            E = in1.efermi - 0.2
        E /= 2.0
        rs[ia] = struct.get_radial(ia)
        vrs[ia] = vsp.get_vsp(ia)
        Zs[ia] = struct.get_z(ia)
        if ia not in los:
            los[ia] = {}
        if l not in los[ia]:
            los[ia][l] = {}
        los[ia][l][idl] = E
    requestial = list(set(requestial))

    # retrive linearization energy of APW basis
    for ia in los:
        apw[ia] = {}
        for l in los[ia]:
            els = in1.get_exceptions(ia)
            E = els[l][0][0]
            if E == 0.3:
                E = in1.efermi - 0.2
            apw[ia][l] = E/2.0

    if opts.debug: 
        # check the linearization energy
        print("LOs:", los)
        print("APW:", apw)

    elos = deepcopy(los)
    for ia, l in requestial:
        # compute APWs
        vr = vrs[ia]
        r = rs[ia]
        Z = Zs[ia]
        apw[ia][l] = list(solve_rde(r, vr, apw[ia][l], l, Z, normalize=True))
        apw[ia][l][0] /= r
        apw[ia][l][1] /= r
        # compute LOs
        idls = list(los[ia][l].keys())
        for idl in idls:
            if idl != 0:
                # LOs
                los[ia][l][idl] = list(solve_rde(r, vr, los[ia][l][idl], l, Z, normalize=True))
            else:
                # lo
                los[ia][l][idl] = list(solve_ene_deriv(r, vr, los[ia][l][idl], l, Z))
            if not opts.uonly:
                print("compute LO/lo coefficients and normalize...")
                coef = - los[ia][l][idl][2] / apw[ia][l][2]
                los[ia][l][idl][0] += coef * apw[ia][l][0] * r
                los[ia][l][idl][1] += coef * apw[ia][l][1] * r
                normsq = logr_int(r, los[ia][l][idl][0], los[ia][l][idl][1], \
                                  los[ia][l][idl][0], los[ia][l][idl][1], )
                los[ia][l][idl][0] /= np.sqrt(normsq)
                los[ia][l][idl][1] /= np.sqrt(normsq)
            sign = np.sign(los[ia][l][idl][0][0])
            los[ia][l][idl][0] /= r * sign
            los[ia][l][idl][1] /= r * sign
    
    # plotting
    for ia, l in requestial:
        r = rs[ia]
        for idl, data in los[ia][l].items():
            typeStr = {0: "lo"}.get(idl, "LO")
            E = elos[ia][l][idl]
            n = data[-1]
            ax.plot(r, data[0], label="Atom=%d l=%d E=%8.4f a.u. (%d nodes) (%s)" % (ia, l, E, n, typeStr), \
                lw=2)
    # adjust Axes attributes
    ax.set_ylabel("u(r)", size=20)
    ax.set_xlabel("r (a.u.)", size=20)
    for a in ['top', 'bottom', 'left', 'right']:
        ax.spines[a].set_linewidth(4)
    ax.axhline(0.0, lw=2, ls=":", color="k")
    ax.tick_params(axis="both", which="both", direction="in", width=2, top=True, right=True)
    ax.tick_params(axis="both", which="major", length=14, labelsize=20)
    ax.tick_params(axis="both", which="minor", length=7)
    ax.set_xlim(0.0)


    plt.legend()
    plt.show()


if __name__ == "__main__":
    pw_plot_apw()