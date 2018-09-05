#!/usr/bin/env python

# ====================================================
#
#     File Name : pw_merge.py
# Creation Date : 20-10-2017
# Last Modified : Fri 20 Oct 2017 03:51:16 PM CST
#    Created By : Min-Ye Zhang
#                 (Based on Xu Xi's script)
#       Contact : stevezhang@pku.edu.cn
#       Purpose : Merge the energy file of wien2k parallel computing
#
# ====================================================

from __future__ import print_function
from w2k_utils import w2k_get
from pw_anal_utils import Get_Casename
import os,sys
import argparse
import subprocess as sp
import fnmatch

def Main(ArgList):

# ==============================

    description='Merge parallel files of wien2k'

    parser=argparse.ArgumentParser(description=description,
    		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f',default='energy',dest='filename',help="the name of files")
    opts=parser.parse_args()

# ==============================

    casename = Get_Casename()
    filename = opts.filename
    
    files_list = fnmatch.fnmatch('./',casename+'.'+filename+'_*')

    if files_list == False:
        nproc = 0
    else:
        nproc = len(files_list)

    
    # skip merging when the calculation is serial
    if nproc==0:
        print("Serial calculation detected. No need to merge.")
    else:
        try:
            print("Attempt to use 'x joinvec' to merge .energy and .vector files...")
            sp.check_output('x joinvec',shell=True)
            return
        except:
            print("Failed to use joinvec... Mannually merge energy files only...")

        nat=w2k_get(casename,'nat') #number of atoms
        sp.check_call('cat %s.%s_1 > %s.%s' %(casename,filename,casename,filename),shell=True)
        energy_file=file('%s.%s' %(casename,filename),'a')
        for i in range(nproc-1):
            append_file=file(casename+'.'+filename+'_'+str(i+2),'r')
            for j in range(nat*2): #skip lines of linearization energy
                append_file.readline()
            for line in append_file:
                energy_file.write(line) #append lines of eigenvalues of kpoints
            append_file.close()
        energy_file.close()

    return

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

