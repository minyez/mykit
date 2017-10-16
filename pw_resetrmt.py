#!/usr/bin/env python
# use setrmt_lapw and clminter to restart the wien2k calculation from a different RMT setting

from pw_anal_utils import Get_Casename
from sys import argv,exit
import subprocess as sp
from shutil import copy2

description = 'Input in the format of ``Element:RMT'' one by one, e.g. Cu:2.1 Cl:2.1.'


case = Get_Casename()
case_struct = case+'.struct'
new_rmt = argv[1:]

if len(argv) == 1:
    print "Current case: " + case
    print description
    exit(1)

for i in xrange(len(new_rmt)):
    new_rmt[i] = new_rmt[i].strip()
    while new_rmt[i].startswith(','):
        new_rmt[i] = new_rmt[i][1:].strip()
    while new_rmt[i].endswith(','):
        new_rmt[i] = new_rmt[i][:-1].strip()

new_rmt_str = ','.join(new_rmt)
sp.call("setrmt_lapw "+case+" -a "+new_rmt_str,shell=True)
copy2(case_struct+'_setrmt',case_struct+'_new')
# for sp case, also need x clminter -up/-dn
sp.call("x clminter",shell=True)

copy2(case_struct+'_new',case_struct)
copy2(case+'.clmsum_new',case+'.clmsum')


