#!/usr/bin/env python

from __future__ import print_function
import os
import commands
import subprocess as sp
from shutil import copy2, rmtree
from pc_utils import common_io_cleandir
from pv_calc_utils import vasp_vasprun_zmy, vasp_io_change_tag
from pv_anal_utils import vasp_anal_get_fund_gap, vasp_anal_get_BM_info

# pylint: disable=undefined-variable
# variables defined from pv_acfdt.py
nproc = TEMP_nproc
vasp_cmd = TEMP_vasp_cmd

for i in range(6):
    i = str(i)
    common_io_cleandir(i)
    if os.path.isfile('INCAR_'+i):
        copy2('INCAR_'+i, i)
    else:
        raise IOError('INCAR_'+i+' not found')
    for x in ["KPOINTS", "POTCAR", "POSCAR"]:
        copy2(x, i)
for i in range(1, 6):
    copy2("KPOINTS_RPA", str(i)+"/KPOINTS")

# Step 0: converge the charge density with large kpoint mesh
os.chdir('0')
copy2('INCAR_0', 'INCAR')
# remove possible HF tags
vasp_io_change_tag('INCAR', 'NKRED')
vasp_io_change_tag('INCAR', 'PRECFOCK')
vasp_io_change_tag('INCAR', 'TIME')
vasp_io_change_tag('INCAR', 'LHFCALC')
vasp_io_change_tag('INCAR', 'HFSCREEN')
vasp_vasprun_zmy(vasp_cmd, 'out.PBE', 'error.PBE')
copy2('INCAR_0', 'INCAR')
vasp_io_change_tag('INCAR', 'ISTART', new_val=1)
vasp_vasprun_zmy(vasp_cmd, 'out', 'error')
copy2('CHGCAR', '../1')
os.chdir('..')

# Step 1: compute the wave functions with fixed charge
os.chdir('1')
os.rename('INCAR_1', 'INCAR')
vasp_vasprun_zmy(vasp_cmd, 'out', 'error')
Band_info = vasp_anal_get_BM_info()
Eg, f_metal = vasp_anal_get_fund_gap(Band_info)
if f_metal: 
    copy2('WAVECAR', '../5')
BANDS_max = int(sp.check_output("awk '/maximum number/ {print $5}' OUTCAR", shell=True))
nbands = BANDS_max - BANDS_max % nproc
copy2('WAVECAR', '../2')
copy2('WAVECAR', '../3')
copy2('CHGCAR', '../3')

os.chdir('../2')
os.rename('INCAR_2', 'INCAR')
vasp_vasprun_zmy(vasp_cmd, 'out', 'error')
# ==================

os.chdir('../3')
os.rename('INCAR_3', 'INCAR')
vasp_io_change_tag('INCAR', 'NBANDS', new_val=nbands)
if f_metal:
    vasp_io_change_tag('INCAR', 'LOPTICS')

vasp_vasprun_zmy(vasp_cmd, 'out', 'error')
copy2('WAVECAR', '../4')
if not f_metal:
    copy2('WAVEDER', '../4')

os.chdir('../4')
os.rename('INCAR_4', 'INCAR')
vasp_io_change_tag('INCAR', 'NBANDS', new_val=nbands)
if not f_metal:
    vasp_io_change_tag('INCAR', 'SIGMA', new_val=Eg/5.0)

vasp_vasprun_zmy(vasp_cmd, 'out', 'error')
if f_metal:
    os.chdir('../5')
    os.rename('INCAR_5', 'INCAR')
    vasp_vasprun_zmy(vasp_cmd, 'out', 'error')
else:
    os.chdir('../')
    rmtree('5/')
