#!/usr/bin/env python

# write wien2k machine file

from argparse import ArgumentParser
from pw_anal_utils import Get_Basename
import subprocess as sp
import os,sys

def read_kgen(kgen_file):
    with open(kgen_file,'r') as kgen:
        try:
            klines = kgen.readlines()
        except IOError:
            print "Should have KGEN file. Exit."
            sys.exit(1)
    nkpt = int(klines[0].split()[0])
    return nkpt

def remove_machine():
    '''
    Remove .machinex files, x integer
    '''
    machine_list = sp.check_output('ls .machine*',shell=True).split()
    for mf in machine_list:
        if mf != ".machines":
            os.remove(mf)

# ================= Parser ===================

parser = ArgumentParser()
parser.add_argument("-n", dest="nproc", help="number of processors to use",type=int, default=1)
parser.add_argument("-D", dest="debug", help="debug mode",action='store_true')
opts = parser.parse_args()

hostname = sp.check_output("hostname",shell=True).split()[0]
if hostname == "MZ":
    hostname = "localhost"
if opts.debug:
    print hostname

dirname = Get_Basename()

if opts.debug:
    print dirname

nproc = opts.nproc
print "Number of processors: ", nproc

# for hybrid functional calculations
if os.path.exists(dirname+'.kgen_ibz'):
    kgen_file = dirname + '.kgen_ibz'
else:
    kgen_file = dirname + '.kgen'

nkpt = read_kgen(kgen_file)
    
if nproc <=1:
    print "Less than 1 processors. No need for parallel calculation."
else:
    if nproc >= nkpt:
        kptlist = nkpt * ['1']
    else:
        k_per_proc = nkpt/nproc
        k_over = nkpt - nproc * k_per_proc
        kptlist = [str(k_per_proc+1)] * k_over + [str(k_per_proc)] * (nproc - k_over)
    if opts.debug:
        print kptlist

print "Writing .machines"
with open('.machines','w') as f:
    for x in kptlist:
        f.write(x+':'+hostname+'\n')
    f.write("\ngranularity:1\nextrafine:1")

remove_machine()
