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

import sys, os
from argparse import ArgumentParser
from pv_calc_utils import vasp_write_incar_minimal_elec,\
                          vasp_io_get_NPAR,vasp_vaspcmd_zmy,\
                          vasp_io_change_tag,vasp_write_incar_exact,\
                          vasp_write_kpoints_basic

# =====================================================

def acfdt_write_incar_nscHF(incar,tag_xc,encut,npar=1):
    incar.write("Non-self-consistent Hatree-Fock\n\n")
    vasp_write_incar_minimal_elec(incar,tag_xc,encut,npar=npar,mode_smear=0,wfrestart=1)
    incar.write("\n# nsc-HF\n")
    incar.write(" AEXX    = 1.0; ALDAC   = 0.0; AGGAC   = 0.0\n")
    incar.write(" ALGO  = EIGENVAL\n")
    incar.write(" LHFCALC = .TRUE.\n")
    incar.write(" LWAVE   = .FALSE.; LCHARG = .FALSE.\n")
    incar.write(" NELM   = 1\n")

# =====================================================

def acfdt_write_incar_rpa(incar,tag_xc,encut,nomega=16,omegatl=1000,mode_smear=[-1,0.01]):
    '''
    Write INCAR for ACFDT-RPA calculation
    '''
    incar.write("ACFDT-RPA\n\n")
    vasp_write_incar_minimal_elec(incar,tag_xc,encut,mode_smear=mode_smear,wfrestart=1)
    incar.write(" NBANDS = 100\n") # will be modified after first step
    incar.write("\n# RPA Frequency integration\n")
    incar.write(" ALGO    = ACFDT\n")
    incar.write(" OMEGATL = %d\n" % omegatl)
    incar.write(" NOMEGA  = %d\n" % nomega)

# =====================================================

def acfdt_write_running_script(vasp_cmd,nproc):
    '''
    Write the python script running RPA
    '''
    module_import = "#!/usr/bin/env python\n" +    \
            "import os, commands\n"  +             \
            "import subprocess as sp\n"  +         \
            "from shutil import copy2,rmtree\n" +  \
            "from pc_utils import common_io_cleandir\n" +  \
            "from pv_calc_utils import vasp_vasprun_zmy,vasp_io_change_tag\n" + \
            "from pv_anal_utils import vasp_anal_get_fund_gap,vasp_anal_get_BM_info\n\n"
    print "Wrinting executing script ..."
    vasp_cmd_list = vasp_cmd.split()
    vasp = vasp_cmd_list[-1]
    if len(vasp_cmd_list) == 0:
        sys.exit("Wrong vasp executives")
    with open("acfdt.py",'w') as ofile:
        ofile.write(module_import)
        ofile.write("nproc = %s\n" % nproc)
# set metagga to True and preconverge PBE orbitals for metaGGA calculation, or instead use CG in the first step
        ofile.write("metagga = False\n")
        ofile.write("vasp = '%s'\n" % vasp)
        ofile.write("vasp_cmd = '%s'\n" % vasp_cmd)
        ofile.write('''for i in xrange(6):
    i = str(i)
    common_io_cleandir(i)
    if os.path.isfile('INCAR_'+i):
        copy2('INCAR_'+i,i)
    for x in ["KPOINTS","POTCAR","POSCAR"]:
        copy2(x,i)\n''')
#        ofile.write("\n# ==================\n\nif \nos.rename('INCAR_1','INCAR')\n")
        ofile.write('''if os.path.isfile('INCAR_0'):
    os.chdir('0')
    copy2('INCAR_0','INCAR')
    vasp_io_change_tag('INCAR','NKRED')
    vasp_io_change_tag('INCAR','PRECFOCK')
    vasp_io_change_tag('INCAR','TIME')
    vasp_io_change_tag('INCAR','LHFCALC')
    vasp_io_change_tag('INCAR','HFSCREEN')
    vasp_vasprun_zmy(vasp_cmd,'out.PBE','error.PBE')
    os.rename('INCAR_0','INCAR')
    vasp_io_change_tag('INCAR','ISTART',new_val=1)
    vasp_vasprun_zmy(vasp_cmd,'out','error')
    copy2('WAVECAR','../1')
    os.chdir('..')''')
        ofile.write("\n# ==================\n\nos.chdir('1')\nos.rename('INCAR_1','INCAR')\n")
        ofile.write("if metagga:vasp_io_change_tag('INCAR','ALGO',new_val = 'A')\n")
# run vasp
        ofile.write("vasp_vasprun_zmy(vasp_cmd,'out','error')\n")
        ofile.write("Band_info = vasp_anal_get_BM_info()\n")
        ofile.write("Eg,f_metal = vasp_anal_get_fund_gap(Band_info)\n")
        ofile.write("if f_metal: copy2('WAVECAR','../5')\n")
        ofile.write('''BANDS_max = int(sp.check_output("awk '/maximum number/ {print $5}' OUTCAR",shell=True))\n''')
        ofile.write("nbands = BANDS_max - BANDS_max % nproc\n")
        ofile.write("copy2('WAVECAR','../2')\n")
        ofile.write("copy2('WAVECAR','../3')\n")
        ofile.write("\n# ==================\n\nos.chdir('../2')\nos.rename('INCAR_2','INCAR')\n")
# run vasp
        ofile.write("vasp_vasprun_zmy(vasp_cmd,'out','error')\n")
        ofile.write("\n# ==================\n\nos.chdir('../3')\nos.rename('INCAR_3','INCAR')\n")
        ofile.write("vasp_io_change_tag('INCAR','NBANDS',new_val = nbands)\n")
        ofile.write("if f_metal: vasp_io_change_tag('INCAR','LOPTICS')\n")
# run vasp
        ofile.write("vasp_vasprun_zmy(vasp_cmd,'out','error')\n")
        ofile.write("copy2('WAVECAR','../4')\n")
        ofile.write("if not f_metal: copy2('WAVEDER','../4')\n")
        ofile.write("\n# ==================\n\nos.chdir('../4')\nos.rename('INCAR_4','INCAR')\n")
        ofile.write("vasp_io_change_tag('INCAR','NBANDS',new_val = nbands)\n")
        ofile.write("if not f_metal: vasp_io_change_tag('INCAR','SIGMA',new_val = Eg/5.0)\n")
# run vasp
        ofile.write("vasp_vasprun_zmy(vasp_cmd,'out','error')\n")
        ofile.write("\n# ==================\n\nif f_metal:\n")
        ofile.write("    os.chdir('../5')\n")
        ofile.write("    os.rename('INCAR_5','INCAR')\n")
# run vasp
        ofile.write("    vasp_vasprun_zmy(vasp_cmd,'out','error')\n")
        ofile.write("else:\n")
        ofile.write("    os.chdir('../')\n")
        ofile.write("    rmtree('5/')")

# =====================================================

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Prepare input files and executing script for ACFDT-RPA
    calculation. Need to first prepare POSCAR and POTCAR
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-e",dest='encut',type=int,default=0,help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors ")
    parser.add_argument("-x",dest='tag_xc',default=None,help="type of XC functional")
    parser.add_argument("-k",dest='nk',type=str,default=[6,6,6],help="Number of kpoints, '[kx,ky,kz]', quotes needed")
    parser.add_argument("-w",dest='nomega',default=16,help="Number of frequency points")
    parser.add_argument("-v",dest='vasp_path',default="vasp",help="Path of vasp executive")

    opts = parser.parse_args()
    npar = vasp_io_get_NPAR(opts.nproc)

# =====================================================
    print opts.nk
    vasp_path, vasp_cmd = vasp_vaspcmd_zmy(np=opts.nproc,vasp_path=opts.vasp_path)
    if not os.path.isfile('KPOINTS'):
        vasp_write_kpoints_basic(opts.nk)

# Step 0: Preconverge the orbital, if hybrid functional is chosen
    if opts.tag_xc == "HSE06" :
        with open('INCAR_0','w') as incar:
            print "Setting INCAR_0: preconvege for hybrid functional"
            incar.write("Standard DFT\n\n")
            vasp_write_incar_minimal_elec(incar,opts.tag_xc,opts.encut,mode_smear=[0,0.05],npar=npar,wfrestart=False)
        vasp_io_change_tag('INCAR_0','PRECFOCK',new_val='Fast')
        # need to set nkred according to the used KPOINTS
        vasp_io_change_tag('INCAR_0','NKRED',new_val=2)


# Step 1: DFT calculation
    with open('INCAR_1','w') as incar:
        print "Setting INCAR_1: std-DFT"
        incar.write("Standard DFT\n\n")
        vasp_write_incar_minimal_elec(incar,opts.tag_xc,opts.encut,mode_smear=[0,0.05],npar=npar,wfrestart=False)
    vasp_io_change_tag('INCAR_1','NELM',new_val=100)
    if opts.tag_xc == 'HSE06':
        vasp_io_change_tag('INCAR_1','ISTART',new_val=1)

# Step 2: HF calculation
    with open('INCAR_2','w') as incar:
        print "Setting INCAR_2: nsc-HF"
        acfdt_write_incar_nscHF(incar,opts.tag_xc,opts.encut,npar)
    if opts.tag_xc == "HSE06":
        vasp_io_change_tag('INCAR_2','ALGO')
        vasp_io_change_tag('INCAR_2','ALGO',new_val='EIGENVAL')
        vasp_io_change_tag('INCAR_2','TIME')
        vasp_io_change_tag('INCAR_2','PRECFOCK')
        vasp_io_change_tag('INCAR_2','HFSCREEN')

# Step 3: Exact diagonalization
#   For RPA calculation of system with a gap, WAVEDER is
#   required to include longwave contribution, while it
#   is better for convergence w.r.t kpoint-mesh if this
#   contribution is neglected.
    with open('INCAR_3','w') as incar:
        print "Setting INCAR_3: Exact Diag."
        vasp_write_incar_exact(incar,opts.tag_xc,opts.encut,npar=npar,mode_smear=[0,0.05],loptics=True)
    if opts.tag_xc == "HSE06":
        vasp_io_change_tag('INCAR_3','ALGO')
        vasp_io_change_tag('INCAR_3','ALGO',new_val='Exact')

# Step 4: ACFDT-RPA calculation
    with open('INCAR_4','w') as incar:
        print "Setting INCAR_4: RPA"
        acfdt_write_incar_rpa(incar,None,opts.encut,nomega=opts.nomega,mode_smear=[-1,0.01])

# Step 5:  HF-correction, if necessary:
    print "Setting INCAR_5: HF-correction (Metal case)"
    vasp_io_change_tag('INCAR_4',"NBANDS",name_ofile="INCAR_5")
    vasp_io_change_tag('INCAR_5','SIGMA',new_val='0.1',backup=False)

# writing executive script
    acfdt_write_running_script(vasp_cmd,opts.nproc)

if __name__ == "__main__":
    Main(sys.argv)

##
