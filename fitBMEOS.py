#!/usr/bin/env python
# Auto non-linear fitting of Birch-Murnaghn Equation of State. 

from __future__ import print_function
import argparse
import numpy as np
from scipy.optimize import curve_fit

# Four argument is needed: data file, number of formula units, E0 and V0. 
# The last two can be extracted by figuring out the lowest point among the data, if not specified.

# pylint: disable=no-member
def BMEOS(V, E0, V0, B0, Bp):
    '''Birch-Murnaghan Equation of State
    '''
    return E0 + 9.0 / 16.0 * V0 * B0 * \
            ((np.power(V0/V, 2.0/3.0) - 1.0)**3*Bp + \
             (np.power(V0/V, 2.0/3.0) - 1.0)**2*(6.0-4.0*np.power(V0/V, 2.0/3.0)))

def readfile(filename, n, startcol=1):
    with open(filename) as f:
        lines = f.readlines()
    i = 0
    vol = []
    ene = []
    while i < len(lines):
        flag = 1
        line = lines[i].split()
        if line[0].startswith("#"):
            i += 1
            continue
        # avoid duplicates
        for v in vol: 
            if v == float(line[startcol].strip()): 
                flag = 0
                break
        if flag:
            vol.append(float(line[startcol].strip())/float(n))
            ene.append(float(line[startcol+1].strip())/float(n))
        i += 1
    data = [np.array(vol), np.array(ene)]
    return data

eVpAcub2GPa = 160.217733

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input data file")
parser.add_argument("-n", help="number of formula units", type=int)
parser.add_argument("-v", help="initial guess of ground state volume", type=float)
parser.add_argument("-e", help="initial guess of ground state energy", type=float)
parser.add_argument("-p", help="flag for plotting", action="store_true")
parser.add_argument("-f", help="fixed B' when fitting (still bug)", type=float)
parser.add_argument("--ncol", type=int, default=3, \
        help="number of data columns. Default is 3: ratio, vol, energy")
args = parser.parse_args()

# initialization
B0 = 1.0
if args.n:
    nfu = args.n
else:
    nfu = 1
if not args.i:
    args.i = "Ener_Vol"
if args.ncol <= 2:
    vol_index = 0
else:
    vol_index = 1


data = readfile(args.i, nfu, startcol=vol_index)

if args.e:
    E0 = args.e
else:
    E0 = min(data[1])
    bottom = np.where(data[1] == E0)
if args.v:
    V0 = args.v
else:
    V0 = data[0][bottom][0]

if not args.f:
    Bp = 4.0
    popt, pcov = curve_fit(BMEOS, data[0], data[1], p0=[E0, V0, B0, Bp])
    Bp = popt[3]
else:
    Bp = args.f
    popt, pcov = curve_fit(lambda V, E0, V0, B0: BMEOS(V, E0, V0, B0, Bp), \
            data[0], data[1], p0=[E0, V0, B0])

print("  E0 = ", popt[0])
print("  V0 = ", popt[1])
print("  B0 = ", popt[2] * eVpAcub2GPa) # modulus, converting eV/A^3 to GPa needs a factor of ~160  
print("  Bp = ", Bp)

# print pcov

# plotting to see if it works
if args.p:
    import pylab
    xfit = np.linspace(0.98*min(data[0]), 1.02*max(data[0]), 100)
    yfit = BMEOS(xfit, popt[0], popt[1], popt[2], Bp)
    pylab.text(0.1, 0.9, 'matplotlib')
    pylab.xlabel("Volume, $\AA^3$")
    pylab.ylabel("Energy, eV")
    pylab.plot(data[0], data[1], marker='o', color='red', ls='')
    pylab.plot(xfit, yfit, color='black')
    pylab.show()

