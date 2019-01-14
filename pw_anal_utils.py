#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pw_anal_utils.py
# Creation Date : 2017-08-24
# Last Modified : Tue 07 Nov 2017 06:59:57 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : provide classes and functions for the analysis
#                 of WIEN2K input and output files
#
# ====================================================

import os,sys
import subprocess as sp
from pc_utils import common_get_dirpath
from fnmatch import fnmatch

# ====================================================

def w2k_get_casename(w2kdir='.'):
    '''
    return the name of the case from case.struct, if there exists a struct file
    or the casename will be set to the name of the directory.
    '''
    w2kdir_abspath = os.path.abspath(w2kdir)
    if os.path.isdir(w2kdir_abspath):
        for filename in os.listdir(w2kdir_abspath):
            if fnmatch(filename, w2kdir_abspath + '/*.struct'):
                case = filename.split('/')[-1][:-8]
                return case

    return os.path.basename(common_get_dirpath(w2kdir))

# ====================================================

def Read_BandStructure(casename,hybrid=False):
    '''
    read the LDA/GGA or HF (hybrid=True) band structure from the case.energy and case.energyhf files, respectively

    :Return
        A list consisting of nkp members, each member is a nband+1-member list, i.e. [[kx,ky,kz], energy of nband]
    '''
#   get case from case.struct
    if casename==None:
        casename = w2k_get_casename()
    print "casename: %s" % casename
    ifile = casename+'.energy'
    Band_Struct = []
    if hybrid or not os.path.exists(ifile):
        print 'Band structure of hybrid functional will be read.'
        ifile = ifile+'hf'
    # exit if case.energyhf is not found
    if not os.path.exists(ifile):
        print '%s not found. Exit.' % ifile
        sys.exit('1')

    with open(ifile,'r') as f:
        lines = f.readlines()

    i = 4
    while i < len(lines):
        words = lines[i].split()
        if len(words) > 2:
        # In case there is no space between kx,ky,kz
            words = lines[i].replace("1-","1 -").replace("0-","0 -").split()
            kp_data = []
            kp_data.append([float(x) for x in words[0:3]])
            nband = int(words[5])
            for j in xrange(i+1,i+nband+1):
                kp_data.append(float(lines[j].split()[1]))
            i += nband+1
            Band_Struct.append(kp_data)
        else:
            i += 1

    return Band_Struct

# ====================================================
