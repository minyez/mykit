#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pv_hf.py
# Creation Date : 10-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function
from pv_calc_utils import vasp_io_get_NPAR, \
                          vasp_vaspcmd_zmy, \
                          vasp_write_incar_minimal_elec, \
                          vasp_io_change_tag
from argparse import ArgumentParser
import os, sys, re
from shutil import copy2

def pv_write_hf_sequence(nproc, vasp_path, incar_pre='INCAR_1', incar_coarse='INCAR_2', incar_stdhf='INCAR_3',seq_out='hf_calc.py'):

    vasp_path, vasp_cmd = vasp_vaspcmd_zmy(np=nproc,vasp_path=vasp_path)
    # check if the three INCARs exists
    incar_list = [incar_pre, incar_coarse, incar_stdhf]
    for i in incar_list:
        if not os.path.exists(i):
            raise ValueError("Error: specified INCAR not found, %s " % i)

    with open(seq_out, 'w') as h_out:
        h_out.write('#!/usr/bin/env python\n')
        h_out.write('# coding=utf-8\n\n')
        h_out.write('from __future__ import print_function\n')
        h_out.write('import os\n')
        h_out.write('from shutil import copy2\n')
        h_out.write('from pv_calc_utils import common_io_cleandir, vasp_vasprun_zmy\n')
        h_out.write('import subprocess as sp\n\n')
        h_out.write('vasp_cmd = "%s"\n\n' % vasp_cmd)
        h_out.write('# check POSCAR, POTCAR and KPOINTS\n')
        h_out.write('ifiles = ["POSCAR","POTCAR","KPOINTS"]\n')
        h_out.write('for ifile in ifiles:\n')
        h_out.write('    if not os.path.exists(ifile):\n')
        h_out.write('       raise IOError("Error: %s not found" % ifile)\n\n')
        h_out.write('# copy input files except INCAR\n')
        h_out.write('workdirs = ["1_PBE_preconv", "2_HF_coarse", "3_HF_normal"]\n')
        h_out.write('for workdir in workdirs:\n')
        h_out.write('    common_io_cleandir(workdir)\n')
        h_out.write('    for ifile in ifiles:\n')
        h_out.write('        copy2(ifile, workdir)\n\n')
        h_out.write('copy2(\'%s\', workdirs[0]+"/INCAR")\n' % incar_pre)
        h_out.write('copy2(\'%s\', workdirs[1]+"/INCAR")\n' % incar_coarse)
        h_out.write('copy2(\'%s\', workdirs[2]+"/INCAR")\n' % incar_stdhf)
        h_out.write('\n# 1 PBE preconverge\n')
        h_out.write('os.chdir(workdirs[0])\n')
        h_out.write('vasp_vasprun_zmy(vasp_cmd, "out","error")\n')
        h_out.write('\n# 2 coarse HF\n')
        h_out.write('os.chdir(\'../\' + workdirs[1])\n')
        h_out.write('try:\n')
        h_out.write('    copy2("../"+workdirs[0]+"/WAVECAR", "WAVECAR")\n')
        h_out.write('except IOError:\n')
        h_out.write('    print("WARNING: WAVECAR in the previous step not found.")\n')
        h_out.write('vasp_vasprun_zmy(vasp_cmd, "out","error")\n')
        h_out.write('\n# 3 standard HF\n')
        h_out.write('os.chdir(\'../\' + workdirs[2])\n')
        h_out.write('try:\n')
        h_out.write('    copy2("../"+workdirs[1]+"/WAVECAR", "WAVECAR")\n')
        h_out.write('except IOError:\n')
        h_out.write('    print("WARNING: WAVECAR in the previous step not found.")\n')
        h_out.write('vasp_vasprun_zmy(vasp_cmd, "out","error")\n')
        h_out.write('os.chdir("../")\n')


def pv_prepare_hf_calculation(ArgList):

    description = '''Prepare the input files and python executive script for a hybrid functional caluclation (PBE0, HSE06, HF).'''

    parser = ArgumentParser(description=description)
    group1 = parser.add_mutually_exclusive_group()

    parser.add_argument("-e",dest='encut',type=int,default=0,help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors ")
    parser.add_argument("-x",dest='tag_xc',default='HSE06',help="Type of hybrid functional. HSE06|PBE0|HF")
    parser.add_argument("-v",dest='vasp_path',default="vasp",help="Path of vasp executive")
    group1.add_argument("--nkred",dest='nkred',type=int,default=1,help="NKRED set in the coarse calculation of HF ")
    group1.add_argument("--nkredxyz",dest='nkredxyz',help="NKREDX/Y/Z set in the coarse calculation of HF, [X,Y,Z]")
    parser.add_argument("--spin",dest='ispin',type=int,default=1,help="Spin-polarization. 1 for nsp and 2 for sp.")

    opts  = parser.parse_args()
    ispin = opts.ispin
    nproc = opts.nproc
    encut = opts.encut

    npar  = vasp_io_get_NPAR(nproc)

    # generate primitive INCAR files
    with open('INCAR_1','w') as incar:
        vasp_write_incar_minimal_elec(incar,'PBE',encut=encut,npar=npar,spin=ispin)
    with open('INCAR_2','w') as incar:
        vasp_write_incar_minimal_elec(incar,opts.tag_xc,encut=encut,npar=npar,spin=ispin)

    copy2('INCAR_2','INCAR_3')

    vasp_io_change_tag('INCAR_2', 'PRECFOCK', new_val='FAST', backup=False)
    vasp_io_change_tag('INCAR_2', 'ISTART', new_val='1', backup=False)
    vasp_io_change_tag('INCAR_3', 'ISTART', new_val='1', backup=False)

    if opts.nkred != 1:
        vasp_io_change_tag('INCAR_2', 'NKRED', new_val=opts.nkred, backup=False)
    if opts.nkredxyz:
        try:
            string = opts.nkredxyz + 'str'
            nkredxyz = [ x.strip() for x in re.split(r'[,\[\]]',opts.nkredxyz)]
            nkredxyz = [ int(x) for x in ' '.join(nkredxyz).split()]

            nkred_opt = ['NKREDX','NKREDY','NKREDZ']
            for i in range(3):
                if nkredxyz[i] !=1 :
                    vasp_io_change_tag('INCAR_2', nkred_opt[i], new_val=nkredxyz[i], backup=False)
        except:
            print(" WARNING: Unsupported input of NKREDX/Y/Z. Pass")

    pv_write_hf_sequence(nproc=opts.nproc, vasp_path=opts.vasp_path)

# =====================================================

if __name__ == "__main__":
    pv_prepare_hf_calculation(sys.argv)

