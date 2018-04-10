#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pc_sync_mykit.py
# Creation Date : 10-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function
import sys, os, re
import subprocess as sp
from argparse import ArgumentParser

def pc_read_scsites(f_scsites=None):
    '''
    Read the path from scsites file.
    Each line in f_scsites should be in the form:
        userid@supercomputerip:path/to/mykit
    '''
    
    # By default, the file storing the supercomputer sites is 'scsites' under the MYKIT directory.
    if f_scsites is None:
        f_scsites = os.path.dirname(__file__)+'/scsites'

    with open(f_scsites,'r') as h_scsites:
        scsites = h_scsites.readlines()

    for i in range(len(scsites)):
        scsites[i] = scsites[i].strip()
        while(scsites[i].endswith('/')):
            scsites[i] = scsites[i][:-1]

    return scsites


def pc_test_connection(scsite):
    '''
    Test if the site is available by pinging it.
    It is not effecitve, since ping is diabled in my local network environment
    '''
    # get the ip address
    site_ip = re.split(r'[@:]', scsite)[1]
    print("Syncing: " + site_ip)


def pc_sync_mykit(ArgList):

    description = '''
    Sync MYKIT to the supercomputer sites, currently by using rsync.
    By default, the file storing the supercomputer sites is 'scsites' under the MYKIT directory.
    You must to set it by yourself and keep it secret.
    PS. for this script, you should upload the SSH public key first to avoid the need of password.
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument('-f',dest='f_scsites',default=None,help="The file storing the supercomputer sites")
    opts = parser.parse_args()

    supercomputer_sites = pc_read_scsites(opts.f_scsites)

    rsync_cmd = "rsync -qazru --inplace "
    local_mykit_path = os.path.dirname(os.path.abspath(__file__))

    #print(supercomputer_sites)
    #print(local_mykit_path)

    for site in supercomputer_sites:
        pc_test_connection(site)
        sp.check_output(rsync_cmd + " " + local_mykit_path + "/ " + site, shell=True)



# ==============================

if __name__ == "__main__":
    pc_sync_mykit(sys.argv)

