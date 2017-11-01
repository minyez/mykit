#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pw_anal_utils.py
# Creation Date : 2017-08-24
# Last Modified : Wed 01 Nov 2017 05:30:11 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : provide classes and functions for the analysis
#                 of WIEN2K input and output files
#
# ====================================================

import os,sys
import subprocess as sp

# ====================================================

def Get_ParentDir(path='PWD'):
    '''
    return the parent directory of path.
    Return: path,                         if path is a directory itself
            the parent directory of path, if path is a file
    '''
    if os.path.isdir(path):
        parentdir = os.path.abspath(path)
    if os.path.isfile(path):
        parentdir = os.path.abspath(path)[:-len(os.path.basename(path))]
    return parentdir

# ==

def Get_Basename():
    '''
    return the name of current directory
    '''
    path = os.environ['PWD']
    dirname = os.path.basename(path)

    return dirname

# ====================================================

def Get_Casename():
    '''
    return the name of the case from case.struct
    '''
    path = sp.check_output('ls *.struct',shell=True)
    case = path[:-8]

    return case

# ====================================================

def Read_BandStructure(hybrid=False):
    '''
    read the LDA/GGA or HF (hybrid=True) band structure from the case.energy and case.energyhf files, respectively

    :Return
        A list consisting of nkp members, each member is a nband+1-member list, i.e. [[kx,ky,kz], energy of nband]
    '''
#   get case from case.struct
    case = Get_Casename()

    ifile = case+'.energy'
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
