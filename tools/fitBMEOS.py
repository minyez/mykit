#!/usr/bin/env python3
'''Auto non-linear fitting of Birch-Murnaghn Equation of State. 

The data file should be prepared forehand, and two formats are supported.

TODO: 
    extend to EOS with number of parameters other than 4
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from mykit.core.constants import EV_PER_ANG_CUB2GPA
from mykit.core.eos import Birch_Murnaghan


def _readfile(filename, v_col, e_col, nfu=1):
    data = np.loadtxt(filename, unpack=True)
    return data[v_col]/nfu, data[e_col]/nfu


def fitBMEOS():
    '''Main stream to fit BM-EOS
    '''
    parser = ArgumentParser(__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", dest="input", type=str, default="Ener_Vol", \
        help="input data file")
    parser.add_argument("-n", dest="nunits", type=int, default=1, \
        help="number of formula units")
    #parser.add_argument("-v", type=float, \
    #    help="initial guess of ground state volume")
    #parser.add_argument("-e", type=float, \
    #    help="initial guess of ground state energy")
    parser.add_argument("-p", dest="plot", action="store_true", \
        help="flag for showing the plot")
    parser.add_argument("--fbp", type=float, \
        help="fixed B' when fitting (still bug)")
    parser.add_argument("-v", dest="v_col", type=int, default=1, \
        help="column of volume")
    parser.add_argument("-e", dest="e_col", type=int, default=2, \
        help="column of energy")
    parser.add_argument("--gpa", dest="gpa", action="store_true", \
        help="Use GPa as unit to print bulk modulus")
    parser.add_argument("--png", dest="png", type=str, default=None, \
        help="Output to png. Filename wo png extension.")
    parser.add_argument("--write-fit", dest="write_fit", action="store_true", \
        help="Write fitted and shifted original data.")
    args = parser.parse_args()

    # initialization
    B0 = 1.0
    nfu = args.nunits

    vs, es = _readfile(args.input, args.v_col-1, args.e_col-1, nfu=nfu)
    opts_curve_fit = {"maxfev": 5000}

    # initial guess
    E0 = np.min(es)
    bottom = np.argmin(es)
    V0 = vs[bottom]
    Bp = 0.0

    if not args.fbp:
        popt, pcov = curve_fit(Birch_Murnaghan, vs, es, p0=[E0, V0, B0, Bp], \
                **opts_curve_fit)
        Bp = popt[3]
    else:
        Bp = args.fbp
        popt, pcov = curve_fit(lambda V, E0, V0, B0: Birch_Murnaghan(V, E0, V0, B0, Bp), \
                vs, es, p0=[E0, V0, B0], **opts_curve_fit)

    perr = np.sqrt(np.diag(pcov))
    gpaConvDict = {True: EV_PER_ANG_CUB2GPA, False: 1.0}
    conv = gpaConvDict[args.gpa]

    print("  E0 = %13.5f (%9.5f)" % (popt[0], perr[0]))
    print("  V0 = %13.4f (%9.4f)" % (popt[1], perr[1]))
    print(" (a0 = %13.4f for FCC)" % np.power(4*popt[1], 1.0/3))
    print("  B0 = %13.3f (%9.3f)" % (popt[2] * conv, perr[2] * conv))
    print("  Bp = ", Bp)

    # compute R2
    SSres = np.sum((es-Birch_Murnaghan(vs, *popt[0:3], Bp))**2)
    SStot = np.sum((es-np.mean(es))**2)
    R2 = 1 - SSres / SStot
    print("  R2 = ", R2)

    xfit = np.linspace(0.98*min(vs), 1.02*max(vs), 100)
    yfit = Birch_Murnaghan(xfit, *popt[0:3], Bp)
    # shift -E0
    yfit_shifted = yfit - popt[0]

    # plotting to see if it works
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.set_xlabel(r"Volume, $\AA^3$")
    ax.set_ylabel("Energy, eV")
    annot = f"R2={R2:.3f}\nE0={popt[0]:.4f}($\pm${perr[0]:.4f})\nV0={popt[1]:.2f}($\pm${perr[1]:.2f})\nB0(GPa)={popt[2]*EV_PER_ANG_CUB2GPA:6.1f}($\pm${perr[2]*EV_PER_ANG_CUB2GPA:.1f})"
    ax.plot(vs, es, marker='o', color='k', ls='none')
    ax.plot(xfit, yfit, color='k')
    ax.annotate(s=annot, xy=(0.20, 0.85), xycoords="axes fraction", fontsize=16)
    fig.tight_layout()
    if args.png is not None:
        plt.savefig(args.png+'.png')
    if args.plot:
        plt.show()

    if args.write_fit:
        # write fitted data
        with open(args.input + "_fitted", 'w') as h:
            h.write("#Fitted Birch-Murnaghan EOS\n")
            h.write("#E0(eV) %s, std deviation %s\n" % (popt[0], perr[0]))
            h.write("#V0(A3) %s, std deviation %s\n" % (popt[1], perr[1]))
            h.write("#B0(eV/A3) %s, std deviation %s\n" % (popt[2], perr[2]))
            h.write("#B0(GPa) %s, std deviation %s\n" % (popt[2] * EV_PER_ANG_CUB2GPA, perr[2] * EV_PER_ANG_CUB2GPA))
            h.write("#B' %s\n" % Bp)
            h.write("#R2 %s\n" % R2)
            for i, fd in enumerate(yfit):
                h.write("%s  %s\n" % (xfit[i], fd))
        # shift the original data 
        with open(args.input + "_shifted", 'w') as h:
            print("#Shifted original data by fitted E0=%s" % popt[0], file=h)
            es -= popt[0]
            for v, e in zip(vs, es):
                print("%s  %s" % (v, e), file=h)
            print("", file=h)
            print("#Shifted fitted data by fitted E0=%s" % popt[0], file=h)
            for i, y in enumerate(yfit_shifted):
                print("%s  %s" % (xfit[i], y), file=h)


if __name__ == '__main__':
    fitBMEOS()
