#!/usr/bin/env python3
'''Auto non-linear fitting of Birch-Murnaghn Equation of State. 

The data file should be prepared forehand, and two formats are supported.

    1. volume ratio, volume, energy (default)
    2. volume, energy

To use the second format, add option '--ncol 2'
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np
from scipy.optimize import curve_fit

from mykit.core.constants import evPerAcub2gpa
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
    parser = ArgumentParser(__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", dest="input", type=str, default="Ener_Vol", \
        help="input data file")
    parser.add_argument("-n", dest="nunits", type=int, default=1, \
        help="number of formula units")
    parser.add_argument("-v", type=float, \
        help="initial guess of ground state volume")
    parser.add_argument("-e", type=float, \
        help="initial guess of ground state energy")
    parser.add_argument("-p", dest="plot", action="store_true", \
        help="flag for plotting")
    parser.add_argument("--fbp", type=float, \
        help="fixed B' when fitting (still bug)")
    parser.add_argument("--ncols", type=int, default=3, \
        help="number of data columns. Default is 3: ratio, vol, energy")
    parser.add_argument("--gpa", dest="gpa", action="store_true", \
        help="Use GPa as unit to print bulk modulus")
    args = parser.parse_args()

    # initialization
    B0 = 1.0
    if args.nunits:
        nfu = args.nunits
    else:
        nfu = 1
    if args.ncols <= 2:
        vol_index = 0
    else:
        vol_index = 1


    data = _readfile(args.input, nfu, startcol=vol_index)

    if args.e:
        E0 = args.e
    else:
        E0 = min(data[1])
        bottom = np.where(data[1] == E0)
    if args.v:
        V0 = args.v
    else:
        V0 = data[0][bottom][0]

    if not args.fbp:
        Bp = 4.0
        popt, _pcov = curve_fit(Birch_Murnaghan, data[0], data[1], p0=[E0, V0, B0, Bp])
        Bp = popt[3]
    else:
        Bp = args.fbp
        popt, _pcov = curve_fit(lambda V, E0, V0, B0: Birch_Murnaghan(V, E0, V0, B0, Bp), \
                data[0], data[1], p0=[E0, V0, B0])

    gpaConvDict = {True: evPerAcub2gpa, False: 1.0}
    conv = gpaConvDict[args.gpa]

    print("  E0 = ", popt[0])
    print("  V0 = ", popt[1])
    print("  B0 = ", popt[2] * conv) # modulus, converting eV/A^3 to GPa needs a factor of ~160  
    print("  Bp = ", Bp)

    # write fitted data
    xfit = np.linspace(0.98*min(data[0]), 1.02*max(data[0]), 100)
    yfit = Birch_Murnaghan(xfit, popt[0], popt[1], popt[2], Bp)
    # shift -E0
    yfit_shifted = yfit - popt[0]

    with open(args.input + "_fitted", 'w') as h:
        h.write("#Fitted Birch-Murnaghan EOS\n")
        h.write("#E0(eV) %s\n" % popt[0])
        h.write("#V0(A3) %s\n" % popt[1])
        h.write("#B0(eV/A3) %s\n" % popt[2])
        h.write("#B' %s\n" % Bp)
        for i, fd in enumerate(yfit):
            h.write("%s  %s  %s\n" % (xfit[i], fd, yfit_shifted[i]))

    # plotting to see if it works
    if args.plot:
        import pylab
        xfit = np.linspace(0.98*min(data[0]), 1.02*max(data[0]), 100)
        yfit = Birch_Murnaghan(xfit, popt[0], popt[1], popt[2], Bp)
        pylab.text(0.1, 0.9, 'matplotlib')
        pylab.xlabel(r"Volume, $\AA^3$")
        pylab.ylabel("Energy, eV")
        pylab.plot(data[0], data[1], marker='o', color='black', ls='')
        pylab.plot(xfit, yfit, color='black')
        pylab.show()
    # shift the original data 
    with open(args.input + "_orig_shifted", 'w') as h:
        print("#Shifted original data by fitted E0", file=h)
        data[1] = data[1] - popt[0]
        for v, e in zip(data[0], data[1]):
            print("%s  %s" % (v, e), file=h)

if __name__ == '__main__':
    fitBMEOS()
