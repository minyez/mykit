#!/usr/bin/env python

import sys, os, re
from shutil import copy2
from pv_calc_utils import *

# =====================================================


# =====================================================

def Main(ArgList):

# =================== Parser ==========================

    description = '''
    Prepare INCAR files and run the GW calculation
    Need to first prepare POSCAR and POTCAR
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument("-e",dest='encut',type=int,default=0,help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors ")
    parser.add_argument("-b",dest='nbands',default=None,help="Number of bands")
    parser.add_argument("-x",dest='tag_xc',default=None,help="type of XC functional for input orbitals, None for LEXCH in POTCAR")
    parser.add_argument("-k",dest='nk',default=[6,6,6],help="Number of kpoints for SCF, [kx,ky,kz]")
    parser.add_argument("--kg",dest='nk_gw',default=[6,6,6],help="Number of kpoints for GW, [kgx,kgy,kgz]")
    parser.add_argument("-m",dest='f_metal',action = 'store_true',help="flag whether system is metal, F-Semicond., T-metal")
    parser.add_argument("-w",dest='nomega',default=50,help="Number of frequency points")
    parser.add_argument("-D",dest='debug',action='store_true',help="Debug mode")
    parser.add_argument("-v",dest='vasp_path',default="vasp",help="Path of vasp executive")
    parser.add_argument("--gw",dest='gw_mode',default="G0W0",help="Self-consistent level of GW")
    parser.add_argument("--restart",dest='chg_start',action='store_true',help="start from calculating CHGCAR")
    parser.add_argument("--wannier",dest='lwannier',action='store_true',help="Use wannier90 to plot band structure")

    opts = parser.parse_args()
    npar = vasp_io_get_NPAR(opts.nproc)
    nks = opts.nk
    nks_gw = opts.nk_gw
    if opts.debug:
        print nks
        print nks_gw
        sys.exit(0)
    if opts.nbands is not None:
        nbands = opts.nbands
    else:
# setting default, i.e. 8 times the number of processors
        nbands = 8*opts.nproc

# grid of density of states
    nedos = 3000
# lmaxmix, change to larger value if f-element present
    lmm = 4
# set vasp variable
    vasp_path, vasp_cmd = vasp_vaspcmd_zmy(opts.nproc,'mpirun',opts.vasp_path)

    chg_dir = '1_charge'
    exact_dir = '2_exact'
    gw_dir = '3_gw'


# 1) run an SCF-calculation for the charge density
    if not os.path.exists(chg_dir) and opts.chg_start:
# if CHGCAR does not exist but want to start from CHGCAR, report error
        print "Charge has not yet been calculated"
        sys.exit(1)

# if CHGCAR calculated and not start from CHGCAR
    elif os.path.exists(chg_dir+'/CHGCAR') and not opts.chg_start:
        print "Charge calculation done. Move on"

# CHGCAR do not exist, or want to redo CHGCAR calculation
    else:
        if not os.path.exists(chg_dir):
            common_io_checkdir(chg_dir)
        else:    
            print "Charge calculation done before. Restart"
            common_io_cleandir(chg_dir)
        copy2('POSCAR',chg_dir)
        copy2('POTCAR',chg_dir)
        os.chdir(chg_dir)

        with open('INCAR','w') as incar:
            vasp_write_incar_minimal_elec(incar,opts.tag_xc,opts.encut,npar=npar,mode_smear=[0,0.05])
        if opts.f_metal:
            vasp_io_change_tag('INCAR','SIGMA',new_val=0.2)
        vasp_io_change_tag('INCAR','LMAXMIX',new_val=lmm)

        vasp_write_kpoints_basic(nks)
        vasp_vasprun_zmy(vasp_cmd,'out','error')
        os.chdir('..')

# 2) run a non-SCF calculation to get desired number of bands
    common_io_cleandir(exact_dir)
    copy2('POSCAR',exact_dir)
    copy2('POTCAR',exact_dir)
    os.chdir(exact_dir)
    vasp_write_kpoints_basic(nks_gw)
    with open('INCAR','w') as incar:
        vasp_write_incar_exact(incar,opts.tag_xc,opts.encut,nb=nbands,npar=npar,mode_smear=[0,0.05],lwannier=opts.lwannier)
    vasp_io_change_tag('INCAR','NEDOS',new_val=nedos)
    if opts.f_metal:
        vasp_io_change_tag('INCAR','SIGMA',new_val=0.2)

# if kpoints for the first and second calculation is the same, then it is not necessary to fix charge density, and the wavefunction in Step 1 can be used
    if not nks == nks_gw:
        copy2('../'+chg_dir+'/CHGCAR','.')
        vasp_io_change_tag('INCAR','ISTART')
        vasp_io_change_tag('INCAR','ICHARG',new_val=11)
        vasp_io_change_tag('INCAR','LMAXMIX',new_val=lmm)
    else:
        copy2('../'+chg_dir+'/WAVECAR','.')
        vasp_io_change_tag('INCAR','ICHARG')
        vasp_io_change_tag('INCAR','ISTART',new_val=1)
        
    vasp_vasprun_zmy(vasp_cmd,'out','error')
    os.chdir('..')

# 3) run a GW calculation with WAVECAR and WAVEDER from Step 2
# need to add more rules
    common_io_cleandir(gw_dir)
    copy2('POSCAR',gw_dir)
    copy2('POTCAR',gw_dir)
    copy2(exact_dir+'/KPOINTS',gw_dir)
    copy2(exact_dir+'/WAVECAR',gw_dir)
    if not opts.f_metal:
        copy2(exact_dir+'/WAVEDER',gw_dir)
    os.chdir(gw_dir)

    with open('INCAR','w') as incar:
        vasp_write_incar_exact(incar,opts.tag_xc,opts.encut,nb=nbands,npar=1,mode_smear=[0,0.05],lwannier=opts.lwannier)
        incar.write(' NEDOS = %s\n' % nedos )
        incar.write(' NOMEGA = %s\n' % opts.nomega )
    if opts.gw_mode == "G0W0":
        vasp_io_change_tag('INCAR','ALGO',new_val='GW0')
    if opts.f_metal:
        vasp_io_change_tag('INCAR','SIGMA',new_val=0.2)

    vasp_vasprun_zmy(vasp_cmd,'out','error')
    os.chdir('..')

# =====================================================

if __name__ == "__main__":
    Main(sys.argv)
