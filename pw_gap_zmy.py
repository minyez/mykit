#!/usr/bin/env python

from __future__ import print_function
from pw_anal_utils import Read_BandStructure, w2k_get_casename
from argparse import ArgumentParser
import sys,os
import subprocess as sp

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Output the fundamental band gap and the direct band gaps at VBM/CBM,
    from case.energy (LDA/GGA) or case.energyhf (Hybrid functional)
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-f",dest='casename',help="casename")
    parser.add_argument("-d",dest='debug',help="flag for debug mode",action='store_true')
    parser.add_argument("--hf",dest='hybrid',help="flag for hybrid functional",action='store_true')

    opts = parser.parse_args()

    Ry2eV = 13.6056917

# =====================================================
    if not opts.casename:
        casename = w2k_get_casename()
    else:
        casename = opts.casename

    Band_Struct = Read_BandStructure(casename,opts.hybrid)

    if opts.debug:
        print(Band_Struct)
        print(len(Band_Struct))

    scf_file = casename+'.scf'
    if os.path.exists(scf_file):
        pass
    else:
        if os.path.exists(scf_file+'2'):
            scf_file =scf_file+'2'
        else:
            print("No SCF files found. Do calculation first.")
            sys.exit(1)
    # get the number of valence electron
    nelec = sp.check_output("grep :NOE %s | tail -1" % scf_file, shell=True).split()[-1]
    try:
        nelec = int(float(nelec))
        print("The number of valence electrons:  ", nelec)
    except ValueError:
        print("Fail to extract NOE (%s)" % nelec)
        sys.exit(1)

#  for now only non-spin-polarized case
    nelec = nelec/2
    nkp = len(Band_Struct)

    VBM = max([Band_Struct[i][nelec] for i in xrange(nkp)])
    VBM_nk = [Band_Struct[i][nelec] for i in xrange(nkp)].index(VBM)
    CBM = min([Band_Struct[i][nelec+1] for i in xrange(nkp)])
    CBM_nk = [Band_Struct[i][nelec+1] for i in xrange(nkp)].index(CBM)

#   Fundamental gap, direct gap at VBM and CBM
#   Original in Ry, converted to eV
    Eg = (CBM -VBM) * Ry2eV
    direct_gap_VBM = (Band_Struct[VBM_nk][nelec+1] - Band_Struct[VBM_nk][nelec]) * Ry2eV
    direct_gap_CBM = (Band_Struct[CBM_nk][nelec+1] - Band_Struct[CBM_nk][nelec]) * Ry2eV

    print(":Band Gap  %6.3f eV " % Eg)
    if VBM_nk == CBM_nk:
        if opts.debug:
            print(Band_Struct[VBM_nk][0])
        print("Direct gap at (%6.3f, %6.3f, %6.3f)" % tuple(Band_Struct[VBM_nk][0]))
    else:
        print("Indirect gap:")
        print(":VBM = (%6.3f, %6.3f, %6.3f) Eg(direct) = %6.3f eV  nk=%i" \
        % (Band_Struct[VBM_nk][0][0], Band_Struct[VBM_nk][0][1], Band_Struct[VBM_nk][0][2], direct_gap_VBM, VBM_nk+1))
        print(":CBM = (%6.3f, %6.3f, %6.3f) Eg(direct) = %6.3f eV  nk=%i" \
        % (Band_Struct[CBM_nk][0][0], Band_Struct[CBM_nk][0][1], Band_Struct[CBM_nk][0][2], direct_gap_CBM, CBM_nk+1))

if __name__ == "__main__":
    Main(sys.argv)
