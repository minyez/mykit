#!/usr/bin/env python
# calculate the DOS up to a particular band index, default is the last band calculated in vasp

from pv_anal_utils import vasp_anal_read_eigen
from argparse import ArgumentParser
import sys,math
import numpy as np
import subprocess as sp


def gaussian(x_array,bande,sigma):
    
    a = 1.0E0/math.sqrt(2.0E0*math.pi)/sigma
    exp = np.exp(-np.power((x_array-bande)/sigma,2)/2.0E0)

    return a * exp


def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Calculate the DOS up to a particular band index. Default is the last band calculated in vasp.
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')
    parser.add_argument("-p",dest='prev',help="flag for preview the DOS",action='store_true')
    parser.add_argument("-g",dest='vbm_zero',help="flag to set the VBM as energy zero. Default the Fermi energy",action='store_true')
    parser.add_argument("-n",dest='bandf',help="the last band till which you want to calculate the DOS",type=int,default=10000)
    parser.add_argument("-s",dest='sigma',help="the Gaussian smearing for DOS (in eV)",type=float,default=0.05)

    opts = parser.parse_args()

# =================== Parameters ======================

#   number of DOS grid points
    ngrid = 3000
#   Gaussian smearing parameter
    sigma = opts.sigma
#   Fermi energy extracted from OUTCAR
    efermi = float(sp.check_output("grep E-fermi OUTCAR | tail -1",shell=True).split()[2])
# =====================================================

    band_struct = vasp_anal_read_eigen(debug=opts.debug)
    if opts.debug:
        print band_struct[0],len(band_struct)

    nelec   = band_struct[0][0]
    nkp     = band_struct[0][1]
    bandmax = band_struct[0][2]

    if opts.bandf == 10000 or opts.bandf > bandmax:
        nbands = bandmax
    else:
        nbands = opts.bandf
    if opts.debug:
        print nbands

    VBM = max([band_struct[ikp+1][nelec/2] for ikp in xrange(nkp)])
    dos_min = min([band_struct[ikp+1][1]      for ikp in xrange(nkp)])
    dos_max = max([band_struct[ikp+1][nbands] for ikp in xrange(nkp)])

    dos_min = dos_min - 5*sigma
    dos_max = dos_max + 5*sigma

    grid_array = np.linspace(dos_min,dos_max,ngrid)
    dos_array = np.zeros(ngrid)

    for ikp in xrange(1,nkp+1):
        for bande in band_struct[ikp][1:nbands+1]:
            dos_array = dos_array + gaussian(grid_array,bande,sigma)

    if opts.vbm_zero:
        grid_array = grid_array - VBM
    else:
        grid_array = grid_array - efermi

    with open('dos.dat','w') as o:
        for i in xrange(ngrid):
            o.write("%8.4f   %15.6f\n" % (grid_array[i],dos_array[i]))
     

    if opts.prev:
        import pylab
        pylab.xlabel("Energy (eV)")
        pylab.ylabel("Density of States (arb. unit)")
        pylab.plot(grid_array,dos_array,color='black')
        pylab.xlim([-4,4])
        pylab.ylim([0,300])
        pylab.show()

# =====================================================

if __name__ == "__main__":
    Main(sys.argv)
