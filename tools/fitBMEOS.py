#!/usr/bin/env python3
'''Auto non-linear fitting of Birch-Murnaghn Equation of State. 

The data file should be prepared forehand, and two formats are supported.

    1. volume, energy (default)
    2. volume ratio, volume, energy

To use the second format, add option '--ncols 3'

TODO: 
    extend to EOS with number of parameters other than 4
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np
from scipy.optimize import curve_fit

from mykit.core.constants import EV_PER_ANG_CUB2GPA
from mykit.core.eos import Birch_Murnaghan


def _readfile(filename, n, startcol=1):
    with open(filename) as f:
        lines = f.readlines()
    _i = 0
    vol = []
    ene = []
    while _i < len(lines):
        flag = 1
        line = lines[_i].split()
        if line[0].startswith("#"):
            _i += 1
            continue
        # avoid duplicates
        for _v in vol: 
            if _v == float(line[startcol].strip()): 
                flag = 0
                break
        if flag:
            vol.append(float(line[startcol].strip())/float(n))
            ene.append(float(line[startcol+1].strip())/float(n))
        _i += 1
    _data = [np.array(vol), np.array(ene)]
    return _data


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
        help="flag for plotting")
    parser.add_argument("--fbp", type=float, \
        help="fixed B' when fitting (still bug)")
    parser.add_argument("--ncols", type=int, default=2, \
        help="number of data columns. Default is 2: vol, energy")
    parser.add_argument("--gpa", dest="gpa", action="store_true", \
        help="Use GPa as unit to print bulk modulus")
    args = parser.parse_args()

    # initialization
    B0 = 1.0
    nfu = args.nunits
    if args.ncols <= 2:
        vol_index = 0
    else:
        vol_index = 1

    data = _readfile(args.input, nfu, startcol=vol_index)
    opts_curve_fit = {"maxfev": 5000}

    # initial guess
    E0 = min(data[1])
    bottom = np.where(data[1] == E0)
    V0 = data[0][bottom][0]
    Bp = 0.0

    if not args.fbp:
        popt, pcov = curve_fit(Birch_Murnaghan, data[0], data[1], p0=[E0, V0, B0, Bp], \
                **opts_curve_fit)
        Bp = popt[3]
    else:
        Bp = args.fbp
        popt, pcov = curve_fit(lambda V, E0, V0, B0: Birch_Murnaghan(V, E0, V0, B0, Bp), \
                data[0], data[1], p0=[E0, V0, B0], **opts_curve_fit)

    perr = np.sqrt(np.diag(pcov))
    gpaConvDict = {True: EV_PER_ANG_CUB2GPA, False: 1.0}
    conv = gpaConvDict[args.gpa]

    print("  E0 = %13.5f (%9.5f)" % (popt[0], perr[0]))
    print("  V0 = %13.4f (%9.4f)" % (popt[1], perr[1]))
    print("  B0 = %13.3f (%9.3f)" % (popt[2] * conv, perr[2] * conv))
    print("  Bp = ", Bp)

    # compute R2
    SSres = np.sum((data[1]-Birch_Murnaghan(data[0], *popt[0:3], Bp))**2)
    SStot = np.sum((data[1]-np.mean(data[1]))**2)
    R2 = 1 - SSres / SStot
    print("  R2 = ", R2)

    xfit = np.linspace(0.98*min(data[0]), 1.02*max(data[0]), 100)
    yfit = Birch_Murnaghan(xfit, *popt[0:3], Bp)
    # shift -E0
    yfit_shifted = yfit - popt[0]

    # plotting to see if it works
    if args.plot:
        import matplotlib.pyplot as plt
        _fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.set_xlabel(r"Volume, $\AA^3$")
        ax.set_ylabel("Energy, eV")
        ax.plot(data[0], data[1], marker='o', color='k', ls='none')
        ax.plot(xfit, yfit, color='k')
        plt.show()

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
            h.write("%s  %s  %s\n" % (xfit[i], fd, yfit_shifted[i]))

    # shift the original data 
    with open(args.input + "_orig_shifted", 'w') as h:
        print("#Shifted original data by fitted E0", file=h)
        data[1] -= popt[0]
        for v, e in zip(data[0], data[1]):
            print("%s  %s" % (v, e), file=h)


if __name__ == '__main__':
    fitBMEOS()
