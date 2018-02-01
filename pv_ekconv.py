#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pv_ekconv.py
# Creation Date : 01-02-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : check the convergence w.r.t kmesh and encut
#                 of a particular system
#
# ====================================================

import sys, os
from shutil import copy2
from datetime import datetime
from argparse import ArgumentParser
from pv_classes import vasp_read_poscar
from pv_anal_utils import vasp_anal_get_enmax
from pv_calc_utils import vasp_vaspcmd_zmy, vasp_vasprun_zmy, vasp_write_incar_minimal_elec, \
                          vasp_io_get_NPAR, vasp_write_kpoints_basic, common_io_checkdir, \
                          vasp_io_change_tag

def Main(ArgList):
    # Parser part
    description = " Check the convergence w.r.t k-point mesh and ENCUT. "
    parser = ArgumentParser(description=description)
    parser.add_argument("-e",dest='encut',type=int,default=0,help="Starting ENCUT. Default 0 will use the largest ENMAX")
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors. Default 1.")
    parser.add_argument("-x",dest='tag_xc',default=None,help="type of XC functional for input orbitals, \
                        None for LEXCH in POTCAR")
    parser.add_argument("--klen",dest='l',type=int,default=30,help="control of kmesh density. Default 30.")
    parser.add_argument("--nec",dest='nencut',type=int,default=6,help="number of encut settings for test. Default 6.")
    parser.add_argument("--nkl",dest='nkleng',type=int,default=6,help="number of klength settings for test. Default 6.")
    parser.add_argument("-m",dest='f_metal',action = 'store_true',help="flag for metal system")
    parser.add_argument("--spin",dest='ispin',type=int, default=1, help="ispin for spin-polarization")
    parser.add_argument("--slab",dest='zdir',type=int, default=0, help="non-periodic direction for slab model. 1|2|3 for a1|a2|a3. Default 0 for 3D periodic.")
    parser.add_argument("-V",dest='view',action='store_true',help="flag for view mode, i.e. view parameters and do not calculate")
    parser.add_argument("-C",dest='clean',action='store_true',help="flag for cleanup of large files, e.g. WAVECAR")
    parser.add_argument("-v",dest='vasp_path',default="vasp",help="path of vasp executive")
    opts = parser.parse_args()

    print " ============ pv_ekconv.py ============"
    # check if POTCAR exists and get the enmax from POTCAR
    if not os.path.exists('POTCAR'):
        print " Error: POTCAR not found. Use pv_addpot.py to generate one. Exit."
        sys.exit(1)
    enmax = vasp_anal_get_enmax()

    # define variables from options
    # Starting ENCUT
    if opts.encut == 0:
        encut_s = enmax
    elif opts.encut < enmax:
        print " Warning: Too small starting ENCUT is set. Use ENMAX instead."
        encut_s = enmax
    else:
        encut_s = opts.encut

    # k length l = k*a
    # ATTENTION: metal would need large l to get converged results
    kleng_s = opts.l

    npar = vasp_io_get_NPAR(opts.nproc)

    # define directory prefix
    encut_dir_prefix = "encut_"
    kleng_dir_prefix = "kleng_"

    # check if INCAR exists. 
    # If not, generate the minimal electronic INCAR according to the options
    if not os.path.exists('INCAR'):
        print " Error: INCAR not found. Generate minimal from options."
        if opts.f_metal:
            mode_smear = 2
        else:
            mode_smear = 0
        with open('INCAR','w') as incar:
            vasp_write_incar_minimal_elec(incar,opts.tag_xc,encut=encut_s,\
                                          mode_smear=mode_smear, npar=npar,spin=opts.ispin)
    if not os.path.exists('POSCAR'):
        print " Error: POSCAR not found. Exit."
        sys.exit(1)

    # generate the ENCUT list and KMESH list
    encut_list = [encut_s + 40.0*i for i in xrange(opts.nencut)]
    kleng_list = [kleng_s +    6*i for i in xrange(opts.nkleng)]
    
    # get vasprun command
    vasp_path, vasp_cmd = vasp_vaspcmd_zmy(opts.nproc, 'mpirun', opts.vasp_path)

    poscar = vasp_read_poscar()
    # get starting time
    dt_s = datetime.now()
    # loop over encut in ENCUT list
    for encut in encut_list:
        print " ---------- ENCUT = %s ----------" % encut
        encut_dir = common_io_checkdir(encut_dir_prefix+str(int(encut)))
        # copy POSCAR POTCAR INCAR
        copy2('POSCAR', encut_dir+'/POSCAR')
        copy2('POTCAR', encut_dir+'/POTCAR')
        copy2('INCAR' , encut_dir+'/INCAR')
        os.chdir(encut_dir)
        # change ENCUT in INCAR to the target value, i.e. encut
        vasp_io_change_tag('INCAR','ENCUT',new_val=encut,backup=False)
        # loop over kleng in klength list
        for kleng in kleng_list:
            # check if this setting is calculated
            kleng_dir = kleng_dir_prefix+str(kleng)
            if os.path.exists(kleng_dir):
                print " - Warning: encut_%s_kleng_%i already calculated. Pass" % (encut,kleng)
                continue
            common_io_checkdir(kleng_dir)
            # copy POSCAR POTCAR INCAR
            copy2('POSCAR', kleng_dir+'/POSCAR')
            copy2('POTCAR', kleng_dir+'/POTCAR')
            copy2('INCAR' , kleng_dir+'/INCAR')
            os.chdir(kleng_dir)
            # write KPOINTS
            nks = [int(kleng/x) for x in poscar.lenlat]
            print " calculating with kmesh: %i %i %i" % (nks[0],nks[1],nks[2])
            vasp_write_kpoints_basic(nks,'G',f_slab=opts.zdir)

            # perform calculation
            if not opts.view:
                vasp_vasprun_zmy(vasp_cmd,'out','error')
            os.chdir('..')
        os.chdir('..')
    # get end time
    dt_e = datetime.now()

    # write executed command in a log file, for future checking
    logfile = 'ekconv.log'
    if not opts.view:
        if os.path.exists(logfile):
            os.rename(logfile,logfile+'_old')
        with open(logfile,'w') as f_log:
            f_log.write("# log file for pv_ekconv.py, executed from %s\n" % dt_s)
            f_log.write("# finished: %s \n" % dt_e)
            f_log.write("# command: %s \n" % " ".join(sys.argv))
            f_log.write("# encut_list: %s \n" % "".join(str(encut_list)))
            f_log.write("# kleng_list: %s \n" % "".join(str(kleng_list)))
        

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

