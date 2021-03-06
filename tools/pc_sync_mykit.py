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
import sys
import os
import re
from argparse import ArgumentParser
import subprocess as sp
from pc_utils import common_print_warn

def common_read_scsites(f_scsites=None):
    '''
    Read the path from scsites file.
    Each line in f_scsites should be in the form:
        userid@supercomputerip:path/to/mykit
    '''
    
    # By default, the file storing the supercomputer sites is 'scsites' under the MYKIT directory.
    if f_scsites is None:
        f_scsites = os.path.dirname(__file__)+'/scsites'

    with open(f_scsites, 'r') as h_scsites:
        scsites = h_scsites.readlines()

    for i, site in enumerate(scsites):
        scsites[i] = scsites[i].strip()
        # remove slash at the end of the site according to the feature of rsync
        while scsites[i].endswith('/'):
            scsites[i] = scsites[i][:-1]

    return scsites


def common_test_connection(scsite):
    '''
    TODO:
        Test if the site is available by pinging it.
        It is not effecitve if ping is diabled in the local network environment
    '''
    # get the ip address
    site_ip = re.split(r'[@:]', scsite)[1]
    print("Syncing: " + site_ip)


def common_sync_once(site, local_mykit_path):

    rsync_cmd = "rsync -qazru --inplace "
    try:
        sp.check_output(rsync_cmd + " " + local_mykit_path + "/p?_*.py " + site, shell=True)
        sp.check_output(rsync_cmd + " " + local_mykit_path + "/README.md " + site, shell=True)
        sp.check_output(rsync_cmd + " " + local_mykit_path + "/templates " + site, shell=True)
        print("Done")
    except sp.CalledProcessError:
        common_print_warn("Error happens in rsync to %s. Pass" % site.split(':')[0], 0)


def common_sync_mykit(ArgList):

    description = '''
    Sync MYKIT to the supercomputer sites, currently by using rsync.
    By default, the file storing the supercomputer sites is 'scsites' under the MYKIT directory.
    You must to set it by yourself and keep it secret.
    PS. for this script, you should upload the SSH public key first to avoid the need of password.
    '''

    parser = ArgumentParser(description=description)
    parser.add_argument('-f', dest='f_scsites', default=None, \
            help="The file storing the supercomputer sites")
    opts = parser.parse_args()

    supercomputer_sites = common_read_scsites(opts.f_scsites)
    local_mykit_path = os.path.dirname(os.path.abspath(__file__))

    #print(supercomputer_sites)
    #print(local_mykit_path)

    for site in supercomputer_sites:
        common_test_connection(site)
        common_sync_once(site, local_mykit_path)



# ==============================

if __name__ == "__main__":
    common_sync_mykit(sys.argv)

