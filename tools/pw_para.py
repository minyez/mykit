#!/usr/bin/env python3

# write wien2k machine file

import subprocess as sp
import os
import socket
from argparse import ArgumentParser
from mykit.wien2k.utils import get_casename

def read_kgen(kgen_file):
    with open(kgen_file, 'r') as kgen:
        try:
            klines = kgen.readlines()
        except FileNotFoundError:
            raise FileNotFoundError("kgen file not found.")
    nkpt = int(klines[0].split()[0])
    return nkpt

def remove_machine():
    '''
    Remove .machinex files, x integer
    '''
    try:
        machine_list = sp.check_output(['ls', '.machine*'], universal_newlines=True).split()
        for mf in machine_list:
            if mf != ".machines":
                os.remove(mf)
    except sp.CalledProcessError:
        pass

def pw_gen_machine_files():

    # ================= Parser ===================
    
    parser = ArgumentParser()
    parser.add_argument("-n", dest="nproc", help="number of processors to use", type=int, default=1)
    parser.add_argument("-D", dest="debug", help="debug mode", action='store_true')
    parser.add_argument("-f", dest="casename", help="case name", default=None)
    opts = parser.parse_args()
    
    #hostname = sp.check_output(["hostname"]).split()[0]
    hostname = socket.gethostname()
    if hostname == "MZ":
        hostname = "localhost"
    if opts.debug:
        print(hostname)
    dirname = opts.casename
    if dirname is None:
        dirname = get_casename()
    if opts.debug:
        print(dirname)
    
    nproc = opts.nproc
    print("Number of processors: ", nproc)
    
    # for hybrid functional calculations
    if os.path.exists(dirname+'.kgen_ibz'):
        kgen_file = dirname + '.kgen_ibz'
    else:
        kgen_file = dirname + '.kgen'
    
    remove_machine()
    nkpt = read_kgen(kgen_file)
    
    if nproc <= 1:
        print("Less than 1 processors. I cannot perform parallel calculation.")
    else:
        if nproc >= nkpt:
            kptlist = nkpt * ['1',]
        else:
            k_per_proc = nkpt//nproc
            k_over = nkpt - nproc * k_per_proc
            kptlist = [str(k_per_proc+1),] * k_over + [str(k_per_proc),] * (nproc - k_over)
        if opts.debug:
            print(kptlist)
    
        print("Writing .machines")
        with open('.machines', 'w') as f:
            for x in kptlist:
                f.write(x+':'+hostname+'\n')
            f.write("\ngranularity:1\nextrafine:1")

if __name__ == "__main__":
    pw_gen_machine_files()
