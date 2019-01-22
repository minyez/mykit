#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pv_acfdt.py
# Creation Date : 2017-04-24
# Last Modified : Fri 01 Dec 2017 07:49:18 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : This script is used to calculate and analyze the ACFDT-RPA results in vasp
#         TO DO :
#                 vasp calculation script : normal dft, HF, exact diagonalization, acfdt
#                                  Analyze: Ener_Vol curve. In another script
#
# ====================================================

from __future__ import print_function
import sys
import os
import re
from shutil import copy2
from argparse import ArgumentParser
from pv_calc_utils import vasp_write_incar_minimal_elec,\
                          vasp_io_get_NPAR, vasp_vaspcmd_zmy,\
                          vasp_io_change_tag, vasp_write_incar_exact,\
                          vasp_write_kpoints_basic

# =====================================================

def acfdt_write_incar_nscHF(incar, tag_xc, encut, npar=1):
    incar.write("Non-self-consistent Hatree-Fock\n\n")
    vasp_write_incar_minimal_elec(incar, tag_xc, encut, npar=npar, \
            mode_smear=0, wfrestart=1)
    incar.write("\n# nsc-HF\n")
    incar.write(" AEXX    = 1.0; ALDAC   = 0.0; AGGAC   = 0.0\n")
    incar.write(" ALGO  = EIGENVAL\n")
    incar.write(" LHFCALC = .TRUE.\n")
    incar.write(" LWAVE   = .FALSE.; LCHARG = .FALSE.\n")
    incar.write(" NELM   = 1\n")

# =====================================================

def acfdt_write_incar_rpa(incar, tag_xc, encut, nomega=16, omegatl=1000, mode_smear=[-1, 0.01]):
    '''
    Write INCAR for ACFDT-RPA calculation
    '''
    incar.write("ACFDT-RPA\n\n")
    vasp_write_incar_minimal_elec(incar, tag_xc, encut, \
            mode_smear=mode_smear, wfrestart=1)
    incar.write(" NBANDS = 100\n") # will be modified after first step
    incar.write("\n# RPA Frequency integration\n")
    incar.write(" ALGO    = ACFDT\n")
    incar.write(" OMEGATL = %d\n" % omegatl)
    incar.write(" NOMEGA  = %d\n" % nomega)

# =====================================================
# pylint: disable=line-too-long
def acfdt_write_running_script(vasp_cmd, nproc):
    '''
    Write the python script running RPA
    '''

    templateName = 'acfdt.py'

    vasp_cmd_list = vasp_cmd.split()
    vasp = vasp_cmd_list[-1]
    if len(vasp_cmd_list) == 0:
        sys.exit("Wrong vasp executives")
    
    cwd = os.environ['PWD']
    copy2(os.path.join(os.path.dirname(__file__), 'templates/'+templateName), cwd)
    os.system('sed -i \'s/TEMP_nproc/%d/g\' %s' % (nproc, templateName))
    os.system('sed -i \'s/TEMP_vasp_cmd/"%s"/g\' %s' % (re.escape(vasp_cmd), templateName))
    

# =====================================================

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Prepare input files and executing script for ACFDT-RPA
    calculation. Need to first prepare POSCAR and POTCAR
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-e", dest='encut', type=int, default=0, \
            help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n", dest='nproc', type=int, default=1, \
            help="Number of processors ")
    parser.add_argument("-x", dest='tag_xc', default=None, \
            help="type of XC functional")
    parser.add_argument("-k", dest='nk', type=str, default=[8, 8, 8], \
            help="Number of kpoints, '[kx,ky,kz]'")
    parser.add_argument("--kgw", dest='nkgw', type=str, default=[4, 4, 4], \
            help="Number of kpoints for RPA, '[kx,ky,kz]'")
    parser.add_argument("-w", dest='nomega', default=16, \
            help="Number of frequency points")
    parser.add_argument("-v", dest='vasp_path', default="vasp", \
            help="Path of vasp executive")

    opts = parser.parse_args()
    npar = vasp_io_get_NPAR(opts.nproc)

# =====================================================
    print(opts.nk)
    vasp_path, vasp_cmd = vasp_vaspcmd_zmy(np=opts.nproc, vasp_path=opts.vasp_path)
    vasp_write_kpoints_basic(opts.nk, koutput="KPOINTS")
    vasp_write_kpoints_basic(opts.nk, koutput="KPOINTS_RPA")

# Step 0: Preconverge the charge 
    with open('INCAR_0', 'w') as incar:
        print("Set INCAR_0 for converging charge density")
        vasp_write_incar_minimal_elec(incar, opts.tag_xc, opts.encut, \
                mode_smear=[0, 0.05], npar=npar, wfrestart=False)
    if opts.tag_xc in ["HSE06", 'PBE0', 'HF', 'SCAN']:
        vasp_io_change_tag('INCAR_0', 'ALGO', new_val='ALL')
    vasp_io_change_tag('INCAR_0', 'LMAXMIX', new_val=6)

# Step 1: DFT calculation for orbitals with fixed charge density
    with open('INCAR_1', 'w') as incar:
        print("Setting INCAR_1: std-DFT")
        incar.write("Standard DFT\n\n")
        vasp_write_incar_minimal_elec(incar, opts.tag_xc, opts.encut, mode_smear=[0, 0.05], \
                npar=npar, wfrestart=False)
    vasp_io_change_tag('INCAR_1', 'NELM', new_val=100)
    vasp_io_change_tag('INCAR_1', 'LORBIT', new_val=11)
    vasp_io_change_tag('INCAR_1', 'LMAXMIX', new_val=6)
    #if opts.tag_xc == 'HSE06':
    #    vasp_io_change_tag('INCAR_1', 'ISTART', new_val=1)

# Step 2: HF calculation
    with open('INCAR_2', 'w') as incar:
        print("Setting INCAR_2: nsc-HF")
        acfdt_write_incar_nscHF(incar, opts.tag_xc, opts.encut, npar)
    if opts.tag_xc == "HSE06":
        vasp_io_change_tag('INCAR_2', 'ALGO')
        vasp_io_change_tag('INCAR_2', 'ALGO', new_val='EIGENVAL')
        vasp_io_change_tag('INCAR_2', 'TIME')
        vasp_io_change_tag('INCAR_2', 'PRECFOCK')
        vasp_io_change_tag('INCAR_2', 'HFSCREEN')

# Step 3: Exact diagonalization
#   For RPA calculation of system with a gap, WAVEDER is
#   required to include longwave contribution, while it
#   is better for convergence w.r.t kpoint-mesh if this
#   contribution is neglected.
    with open('INCAR_3', 'w') as incar:
        print("Setting INCAR_3: Exact Diag.")
        vasp_write_incar_exact(incar, opts.tag_xc, opts.encut, npar=npar, mode_smear=[0, 0.05], loptics=True)
    vasp_io_change_tag('INCAR_3', 'LORBIT', new_val=11)
    vasp_io_change_tag('INCAR_3', 'LMAXMIX', new_val=6)
    if opts.tag_xc == "HSE06":
        vasp_io_change_tag('INCAR_3', 'ALGO')
        vasp_io_change_tag('INCAR_3', 'ALGO', new_val='Exact')

# Step 4: ACFDT-RPA calculation
    with open('INCAR_4', 'w') as incar:
        print("Setting INCAR_4: RPA")
        acfdt_write_incar_rpa(incar, None, opts.encut, nomega=opts.nomega, mode_smear=[-1, 0.01])

# Step 5:  HF-correction, if necessary:
    print("Setting INCAR_5: HF-correction (Metal case)")
    vasp_io_change_tag('INCAR_4', "NBANDS", name_ofile="INCAR_5")
    vasp_io_change_tag('INCAR_5', 'SIGMA', new_val='0.1', backup=False)

# writing executive script
    acfdt_write_running_script(vasp_cmd, opts.nproc)

if __name__ == "__main__":
    Main(sys.argv)

##
