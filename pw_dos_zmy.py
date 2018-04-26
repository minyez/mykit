#!/usr/bin/env python

from pw_anal_utils import Read_BandStructure, w2k_get_casename
from argparse import ArgumentParser
import sys,os
import subprocess as sp

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Ouput the density of states (DOS) by reading the eigenvalues
    from case.energy (LDA/GGA) or case.energyhf (Hybrid functional)
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-d",dest='dosir',help="the directory for saving dos.dat",default='dos')
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')
    parser.add_argument("--hf",dest='hybrid',help="flag for hybrid functional",action='store_true')
    parser.add_argument("-s",dest='smear',help="Parameter for Gaussian smearing",default='0.011')
    parser.add_argument("-g",dest='vbm_zero',help="flag to set the vbm as energy zero",action='store_true')

    opts = parser.parse_args()

    Ry2eV = 13.6056917

# =====================================================

    case = w2k_get_casename()
    Band_Struct = Read_BandStructure(opts.hybrid)
    if opts.debug:
        print Band_Struct
        print len(Band_Struct)

    scf_file = case+'.scf'
    if os.path.exists(scf_file):
        pass
    else:
        if os.path.exists(scf_file+'2'):
            scf_file =scf_file+'2'
        else:
            print "No SCF files found. Do calculation first."
            sys.exit(1)
    efermi = sp.check_output("grep :FER %s | tail -1" % scf_file, shell=True).split()[-1]
    nelec = sp.check_output("grep :NOE %s | tail -1" % scf_file, shell=True).split()[-1]
    try:
        efermi = float(efermi)
        print "E-fermi:  ", efermi
    except ValueError:
        print "Fail to extract E-fermi (%s)" % efermi
        sys.exit(1)
    try:
        nelec = int(float(nelec))
    except ValueError:
        print "Fail to extract NELEC (%s)" % nelec
        sys.exit(1)

#  for now only non-spin-polarized case
    nelec = nelec/2
    nkp = len(Band_Struct)
    VBM = max([Band_Struct[i][nelec] for i in xrange(nkp)])
    dos_min = min([Band_struct[i][0] for i in xrange(nkp)])
    dos_max = max([Band_struct[i][-1] for i in xrange(nkp)])



if __name__ == "__main__":
    Main(sys.argv)
